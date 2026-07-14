// ============================================================================
//  toTrees.cpp : Export a network-from-spr PanMAN as a tskit tree sequence.
//
//  Only the `network-from-spr` architecture is supported: a single base tree
//  (TreeGroup.trees[0]) whose per-block topology is modified by recombination
//  `networkEdges`. Each block is one window's MSA; its consensus is a GAPPED
//  alignment of length `blockLength` (gap = code 0), the block is turned on for
//  every sample by a root block-insertion, and per-node substitutions are
//  stored at MSA columns (nucGapPosition == -1). See src/network.cpp
//  (runFitchForBlock) for how these blocks are produced.
//
//  Coordinates are UNGAPPED positions along the reference tip (default GRCh38),
//  concatenated across the blocks of each chromosome. The reference row is
//  reconstructed along its lineage (consensus + the reference's own nucleotide
//  mutations, including indels) so that its gap structure - and therefore the
//  coordinate map and the reported sequence length - matches the reference's
//  own ungapped length.
//
//  By default only the tables required to build a .trees are written:
//  nodes, edges, and sequence_length. The sites/mutations tables (point
//  substitutions NS / NSNPS at reference-base columns) are large and optional;
//  pass emitMutations=true (CLI: --with-mutations) to also write them. Block
//  mutations and nucleotide indels are never emitted as mutations.
//
//  The output is tskit's documented text-table format. Assemble the final
//  .trees file with scripts/panman_to_trees.py.
// ============================================================================

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panmanUtils {

namespace {

// A nucleotide edit on a node, resolved to an MSA column within its block.
// `isSnv` marks base->base substitutions (NS/NSNPS); indels have isSnv=false
// and are only used when reconstructing the reference row.
struct SubRec {
    int32_t nodeIdx;
    int32_t col;      // MSA column within the block (nucGapPosition assumed -1)
    char    ch;       // resulting character ('-' for deletions)
    bool    isSnv;
};

// A block mutation on a node (used to resolve reference block existence).
struct BMutRec {
    int32_t nodeIdx;
    bool    isBI;
    bool    inversion;
};

// Decode `length` MSA columns of a packed (gapped) block consensus. Unlike
// PanMAN's ungapped reconstruction, code 0 is treated as a gap character here,
// which is how network-from-spr stores its alignment columns.
std::string decodeGappedConsensus(const panmanUtils::Block& block, int64_t length,
                                  panmanUtils::Alphabet alphabet) {
    std::string out;
    if (length <= 0) return out;
    out.reserve(static_cast<size_t>(length));
    const uint8_t perWord = panmanUtils::consensusSymbolsPerPackedWord(alphabet);
    for (int64_t c = 0; c < length; c++) {
        const size_t word = static_cast<size_t>(c) / perWord;
        const size_t idx = static_cast<size_t>(c) % perWord;
        if (word >= block.consensusSeq.size()) { out.push_back('-'); continue; }
        const int code = static_cast<int>(
            panmanUtils::packedConsensusSymbolAt(block.consensusSeq[word], idx, alphabet));
        out.push_back(panmanUtils::getSymbolFromCode(code, alphabet));
    }
    return out;
}

// Extract nucleotide edits from a NucMut. Only main-string edits
// (nucGapPosition == -1) are kept, which is always the case for SPR blocks.
void collectSubs(const panmanUtils::NucMut& m, panmanUtils::Alphabet alphabet,
                 int32_t nodeIdx, std::vector<SubRec>& out) {
    if (m.nucGapPosition != -1) return;
    const uint32_t type = (m.mutInfo & 0x7);
    const int len = (m.mutInfo >> 4);
    switch (type) {
    case panmanUtils::NucMutationType::NS:      // base -> base run
        for (int j = 0; j < len; j++) {
            char ch = panmanUtils::getSymbolFromCode(
                static_cast<int>(panmanUtils::packedMutationSymbolAt(m.nucs, j, alphabet)),
                alphabet);
            out.push_back({ nodeIdx, m.nucPosition + j, ch, true });
        }
        break;
    case panmanUtils::NucMutationType::NI:      // gap -> base run (insertion)
        for (int j = 0; j < len; j++) {
            char ch = panmanUtils::getSymbolFromCode(
                static_cast<int>(panmanUtils::packedMutationSymbolAt(m.nucs, j, alphabet)),
                alphabet);
            out.push_back({ nodeIdx, m.nucPosition + j, ch, false });
        }
        break;
    case panmanUtils::NucMutationType::ND:      // base -> gap run (deletion)
        for (int j = 0; j < len; j++) out.push_back({ nodeIdx, m.nucPosition + j, '-', false });
        break;
    case panmanUtils::NucMutationType::NSNPS:   // single base -> base
        out.push_back({ nodeIdx, m.nucPosition,
                        panmanUtils::getSymbolFromCode(
                            static_cast<int>(panmanUtils::singleMutationSymbol(m.nucs, alphabet)),
                            alphabet),
                        true });
        break;
    case panmanUtils::NucMutationType::NSNPI:   // single gap -> base
        out.push_back({ nodeIdx, m.nucPosition,
                        panmanUtils::getSymbolFromCode(
                            static_cast<int>(panmanUtils::singleMutationSymbol(m.nucs, alphabet)),
                            alphabet),
                        false });
        break;
    case panmanUtils::NucMutationType::NSNPD:   // single base -> gap
        out.push_back({ nodeIdx, m.nucPosition, '-', false });
        break;
    default:
        break;
    }
}

// Pick the reference tip: explicit name, else "GRCh38", else first leaf.
std::string resolveReference(const panmanUtils::Tree& tree,
                             const std::string& requested) {
    if (!requested.empty()) {
        auto it = tree.allNodes.find(requested);
        if (it != tree.allNodes.end()) return requested;
        std::cerr << "Warning: requested reference '" << requested
                  << "' not found; falling back." << std::endl;
    }
    if (tree.allNodes.find("GRCh38") != tree.allNodes.end()) return "GRCh38";
    for (const auto& kv : tree.allNodes) {
        if (kv.second->children.empty()) return kv.first;
    }
    return tree.root ? tree.root->identifier : std::string();
}

// Longest path to a leaf over the union of base-tree child links and network
// child links. Used as node "time" so that time[parent] > time[child] for
// every edge that will be emitted.
std::unordered_map<std::string, int> computeNodeTimes(
    const panmanUtils::Tree& tree,
    const std::vector<panmanUtils::SerializedNetworkEdge>& networkEdges) {

    std::unordered_map<std::string, std::vector<std::string>> children;
    children.reserve(tree.allNodes.size());
    for (const auto& kv : tree.allNodes) {
        const panmanUtils::Node* n = kv.second;
        std::vector<std::string> ch;
        ch.reserve(n->children.size());
        for (const panmanUtils::Node* c : n->children) ch.push_back(c->identifier);
        children[kv.first] = std::move(ch);
    }
    for (const auto& e : networkEdges) {
        if (tree.allNodes.count(e.parentNodeId) && tree.allNodes.count(e.childNodeId)) {
            children[e.parentNodeId].push_back(e.childNodeId);
        }
    }

    std::unordered_map<std::string, int> height;
    std::unordered_map<std::string, int> color; // 0 unvisited, 1 active, 2 done
    height.reserve(tree.allNodes.size());
    color.reserve(tree.allNodes.size());
    bool warnedCycle = false;

    for (const auto& kv : tree.allNodes) {
        if (color[kv.first] == 2) continue;
        std::vector<std::string> stack{ kv.first };
        while (!stack.empty()) {
            const std::string cur = stack.back();
            int& col = color[cur];
            if (col == 0) {
                col = 1;
                for (const std::string& c : children[cur]) {
                    if (color[c] == 0) {
                        stack.push_back(c);
                    } else if (color[c] == 1 && !warnedCycle) {
                        warnedCycle = true;
                        std::cerr << "Warning: cycle detected while assigning node times; "
                                     "some network edges may violate tskit time ordering."
                                  << std::endl;
                    }
                }
            } else if (col == 1) {
                int h = 0;
                for (const std::string& c : children[cur]) {
                    if (color[c] == 2) h = std::max(h, height[c] + 1);
                }
                height[cur] = h;
                col = 2;
                stack.pop_back();
            } else {
                stack.pop_back();
            }
        }
    }
    return height;
}

// Writer for one tree-sequence (a set of four text tables). Nodes are shared
// across all writers; sites/edges/mutations accumulate per writer.
struct TableWriter {
    std::string prefix;          // full path prefix, e.g. ./info/foo.chr1
    std::map<int64_t, char> sites;                       // position -> ancestral
    std::vector<std::tuple<int64_t, int, char>> muts;    // (position, nodeId, derived)
    std::vector<std::tuple<int64_t, int64_t, int, int>> edges; // (left,right,parent,child)

    // Per-child open edge interval for run-length merging across blocks.
    struct OpenEdge { int parent; int64_t left; int64_t right; bool open = false; };
    std::unordered_map<int, OpenEdge> openEdges;

    void addEdge(int child, int parent, int64_t left, int64_t right) {
        OpenEdge& oe = openEdges[child];
        if (oe.open && oe.parent == parent && oe.right == left) {
            oe.right = right;
            return;
        }
        if (oe.open) edges.emplace_back(oe.left, oe.right, oe.parent, child);
        oe.parent = parent;
        oe.left = left;
        oe.right = right;
        oe.open = true;
    }

    void flushEdges() {
        for (auto& kv : openEdges) {
            OpenEdge& oe = kv.second;
            if (oe.open) edges.emplace_back(oe.left, oe.right, oe.parent, kv.first);
            oe.open = false;
        }
    }

    void addSite(int64_t pos, char ancestral) {
        auto it = sites.find(pos);
        if (it == sites.end()) sites[pos] = ancestral;
    }

    void addMut(int64_t pos, int nodeId, char derived) {
        muts.emplace_back(pos, nodeId, derived);
    }
};

inline bool isBase(char c) { return c != '-' && c != 'x' && c != 'N'; }

}  // namespace

void writeTreeSequenceTables(panmanUtils::TreeGroup& TG, const std::string& outPrefix,
                             const std::string& referenceName, bool perChromosome,
                             bool emitMutations) {
    if (TG.trees.empty()) {
        std::cerr << "Error: no trees in PanMAN; cannot export tree sequence." << std::endl;
        return;
    }
    if (TG.trees.size() != 1) {
        std::cerr << "Warning: --to-trees supports the network-from-spr architecture "
                     "(a single base tree). Using trees[0] and ignoring the rest."
                  << std::endl;
    }
    panmanUtils::Tree& tree = TG.trees[0];
    if (!tree.root) {
        std::cerr << "Error: base tree has no root." << std::endl;
        return;
    }
    const panmanUtils::Alphabet alphabet = tree.alphabet;
    const std::string rootId = tree.root->identifier;
    const bool tsDebug = (std::getenv("PANMAN_TS_DEBUG") != nullptr);
    const bool diagOnly = (std::getenv("PANMAN_TS_DIAG") != nullptr);

    const std::string refName = resolveReference(tree, referenceName);
    if (refName.empty() || !tree.allNodes.count(refName)) {
        std::cerr << "Error: could not resolve a reference tip." << std::endl;
        return;
    }
    std::cout << "Using reference tip: " << refName << std::endl;

    // ----- Node ids and times -----
    std::vector<std::string> idToNode;
    std::unordered_map<std::string, int> nodeToId;
    idToNode.reserve(tree.allNodes.size());
    nodeToId.reserve(tree.allNodes.size());
    for (const auto& kv : tree.allNodes) {
        nodeToId[kv.first] = static_cast<int>(idToNode.size());
        idToNode.push_back(kv.first);
    }
    std::unordered_map<std::string, int> nodeTime =
        computeNodeTimes(tree, TG.networkEdges);

    // ----- Reference lineage: nodeIdx -> rank (0 = root ... deepest = reference)
    // The reference is assumed non-recombinant, so its lineage follows the base
    // tree (network overrides are not applied to the reference row itself).
    std::unordered_map<int, int> refPathRank;
    {
        std::vector<int> up;  // reference .. root
        const panmanUtils::Node* n = tree.allNodes.at(refName);
        while (n) {
            up.push_back(nodeToId[n->identifier]);
            if (n->identifier == rootId) break;
            n = n->parent;
        }
        for (size_t i = 0; i < up.size(); i++)
            refPathRank[up[i]] = static_cast<int>(up.size() - 1 - i);
    }

    // ----- Block index (main blocks only) -----
    std::unordered_map<int, const panmanUtils::Block*> blockById;
    blockById.reserve(tree.blocks.size());
    for (const auto& b : tree.blocks) {
        if (b.secondaryBlockId == -1) blockById[b.primaryBlockId] = &b;
    }

    if (tsDebug) {
        std::cerr << "[debug] blockById size=" << blockById.size()
                  << " tree.blocks=" << tree.blocks.size()
                  << " chrList=" << tree.chrList.size() << std::endl;
        for (const auto& chr : tree.chrList) {
            std::cerr << "[debug] chr '" << chr.chrName << "' rawBlockIds:";
            for (size_t i = 0; i < chr.blockIds.size() && i < 12; i++)
                std::cerr << " " << chr.blockIds[i] << "(>>32=" << (chr.blockIds[i] >> 32) << ")";
            std::cerr << std::endl;
        }
    }

    // ----- Chromosome -> ordered primary-block ids -----
    // chrList block ids are normally the encoded (primaryBlockId << 32); some
    // files store them unshifted. Detect which form yields distinct, valid
    // block ids (shifted values collapse to 0 when stored unshifted).
    bool useShift = true;
    if (!tree.chrList.empty()) {
        std::set<int> seen;
        for (const auto& chr : tree.chrList) {
            for (int64_t enc : chr.blockIds) {
                int v = static_cast<int>(enc >> 32);
                if (!blockById.count(v) || !seen.insert(v).second) { useShift = false; break; }
            }
            if (!useShift) break;
        }
    }
    auto decodeBlockId = [&](int64_t enc) -> int {
        return useShift ? static_cast<int>(enc >> 32)
                        : static_cast<int>(enc & 0xFFFFFFFF);
    };

    std::vector<std::pair<std::string, std::vector<int>>> chromBlocks;
    if (!tree.chrList.empty()) {
        for (const auto& chr : tree.chrList) {
            std::vector<int> blocks;
            blocks.reserve(chr.blockIds.size());
            for (int64_t enc : chr.blockIds) blocks.push_back(decodeBlockId(enc));
            chromBlocks.emplace_back(chr.chrName.empty() ? ("chr" + std::to_string(chr.chrIdx))
                                                         : chr.chrName,
                                     std::move(blocks));
        }
    } else {
        std::vector<int> blocks;
        for (const auto& b : tree.blocks)
            if (b.secondaryBlockId == -1) blocks.push_back(b.primaryBlockId);
        std::sort(blocks.begin(), blocks.end());
        chromBlocks.emplace_back("chr0", std::move(blocks));
    }
    if (tsDebug) std::cerr << "[debug] useShift=" << useShift << std::endl;

    // ----- Network-edge overrides per block: block -> (childId -> parentId) -----
    std::unordered_map<int, std::unordered_map<std::string, std::string>> overrides;
    for (const auto& e : TG.networkEdges) {
        for (int64_t bid : e.activeBlockIds)
            overrides[static_cast<int>(bid)][e.childNodeId] = e.parentNodeId;
    }

    // ----- Substitutions and block mutations grouped by block -----
    // To bound memory we keep (a) every SNV only when the mutations table is
    // requested (for sites/mutations) and (b) all edits on the reference lineage
    // (always needed to reconstruct the reference row / coordinate map).
    std::unordered_map<int, std::vector<SubRec>> blockSubs;
    std::unordered_map<int, std::vector<BMutRec>> blockBMuts;
    {
        std::vector<SubRec> tmp;
        for (const auto& kv : tree.allNodes) {
            const int nodeIdx = nodeToId[kv.first];
            const bool onPath = refPathRank.count(nodeIdx) != 0;
            for (const auto& m : kv.second->nucMutation) {
                tmp.clear();
                collectSubs(m, alphabet, nodeIdx, tmp);
                if (tmp.empty()) continue;
                auto& vec = blockSubs[m.primaryBlockId];
                for (const SubRec& s : tmp)
                    if ((emitMutations && s.isSnv) || onPath) vec.push_back(s);
            }
            if (onPath) {
                for (const auto& bm : kv.second->blockMutation) {
                    if (bm.secondaryBlockId != -1) continue;
                    blockBMuts[bm.primaryBlockId].push_back(
                        { nodeIdx, bm.blockMutInfo, bm.inversion });
                }
            }
        }
    }

    // Reconstruct the reference row for one block (gapped consensus + the
    // reference lineage's own edits) and map each non-gap column to an ungapped
    // coordinate. Returns the number of reference bases; fills colCoord (global
    // coord, -1 where the reference is a gap) and refChar (reconstructed row).
    auto reconstructRefBlock = [&](int blockId, int64_t startOffset,
                                   std::vector<int64_t>& colCoord,
                                   std::string& refChar) -> int64_t {
        auto bIt = blockById.find(blockId);
        if (bIt == blockById.end()) return 0;
        const panmanUtils::Block* blk = bIt->second;
        const int64_t L = blk->blockLength;
        if (L <= 0) return 0;

        // Block existence along the reference lineage (root -> reference).
        bool exists = false;
        auto bmIt = blockBMuts.find(blockId);
        if (bmIt != blockBMuts.end()) {
            std::vector<const BMutRec*> onPath;
            for (const BMutRec& bm : bmIt->second)
                if (refPathRank.count(bm.nodeIdx)) onPath.push_back(&bm);
            std::sort(onPath.begin(), onPath.end(), [&](const BMutRec* a, const BMutRec* b) {
                return refPathRank[a->nodeIdx] < refPathRank[b->nodeIdx];
            });
            for (const BMutRec* bm : onPath) {
                if (bm->isBI) exists = true;
                else if (!bm->inversion) exists = false;
            }
        }
        if (!exists) return 0;

        refChar = decodeGappedConsensus(*blk, L, alphabet);

        // Apply the reference lineage's own nucleotide edits (root -> reference).
        auto sIt = blockSubs.find(blockId);
        if (sIt != blockSubs.end()) {
            std::vector<const SubRec*> onPath;
            for (const SubRec& s : sIt->second)
                if (refPathRank.count(s.nodeIdx)) onPath.push_back(&s);
            std::sort(onPath.begin(), onPath.end(), [&](const SubRec* a, const SubRec* b) {
                return refPathRank[a->nodeIdx] < refPathRank[b->nodeIdx];
            });
            for (const SubRec* s : onPath)
                if (s->col >= 0 && s->col < L) refChar[s->col] = s->ch;
        }

        colCoord.assign(static_cast<size_t>(L), -1);
        int64_t running = startOffset;
        for (int64_t c = 0; c < L; c++)
            if (isBase(refChar[c])) colCoord[c] = running++;
        return running - startOffset;
    };

    // ----- Emit -----
    std::filesystem::create_directory("./info");

    // A block id may appear more than once in chrList (observed in some files);
    // each block must contribute its interval only once.
    std::unordered_set<int> processedBlocks;

    auto processChromSet = [&](TableWriter& w,
                               const std::vector<std::pair<std::string, std::vector<int>>>& sets,
                               int64_t& genomeLen) {
        int64_t offset = 0;
        int64_t dbgSeen = 0, dbgExists = 0;
        std::vector<int64_t> colCoord;
        std::string refChar;
        std::string dbgPrevTail;
        int dbgPrevBlock = -1;
        for (const auto& cs : sets) {
            dbgPrevTail.clear();
            dbgPrevBlock = -1;
            for (int blockId : cs.second) {
                if (!blockById.count(blockId)) continue;
                if (!processedBlocks.insert(blockId).second) {
                    if (tsDebug)
                        std::cerr << "[debug]   skipping duplicate block=" << blockId << std::endl;
                    continue;
                }
                dbgSeen++;

                colCoord.clear();
                refChar.clear();
                const int64_t blockRefLen =
                    reconstructRefBlock(blockId, offset, colCoord, refChar);
                if (tsDebug) {
                    const panmanUtils::Block* blk = blockById[blockId];
                    std::string ung;
                    ung.reserve(refChar.size());
                    for (char c : refChar) if (isBase(c)) ung.push_back(c);
                    int ov = 0;
                    if (dbgPrevBlock >= 0 && !dbgPrevTail.empty() && !ung.empty()) {
                        const int cap = std::min<int>(dbgPrevTail.size(), ung.size());
                        for (int k = cap; k >= 50; k--) {
                            if (dbgPrevTail.compare(dbgPrevTail.size() - k, k, ung, 0, k) == 0) { ov = k; break; }
                        }
                    }
                    std::cerr << "[debug]   block=" << blockId
                              << " blockLength=" << blk->blockLength
                              << " refUngapped=" << ung.size()
                              << " overlapWithPrev>=" << ov << std::endl;
                    dbgPrevTail = ung.size() > 200000 ? ung.substr(ung.size() - 200000) : ung;
                    dbgPrevBlock = blockId;
                }
                if (blockRefLen == 0) continue;
                dbgExists++;
                if (diagOnly) { offset += blockRefLen; continue; }
                const int64_t left = offset;
                const int64_t right = offset + blockRefLen;
                offset = right;

                // Edges: marginal parent for each non-root node over [left,right).
                const auto ovIt = overrides.find(blockId);
                for (const auto& kv : tree.allNodes) {
                    if (kv.first == rootId) continue;
                    std::string parent;
                    if (ovIt != overrides.end()) {
                        auto o = ovIt->second.find(kv.first);
                        if (o != ovIt->second.end()) parent = o->second;
                    }
                    if (parent.empty()) {
                        if (!kv.second->parent) continue;
                        parent = kv.second->parent->identifier;
                    }
                    w.addEdge(nodeToId[kv.first], nodeToId[parent], left, right);
                }

                // Sites + mutations: point substitutions at reference columns.
                // Skipped entirely unless the mutations table was requested.
                auto sIt = emitMutations ? blockSubs.find(blockId) : blockSubs.end();
                if (sIt != blockSubs.end()) {
                    for (const SubRec& s : sIt->second) {
                        if (!s.isSnv) continue;
                        if (s.col < 0 || s.col >= static_cast<int32_t>(colCoord.size())) continue;
                        const int64_t coord = colCoord[s.col];
                        if (coord < 0) continue;   // reference is a gap here
                        char anc = refChar[s.col];
                        if (!isBase(anc)) anc = 'N';
                        w.addSite(coord, anc);
                        w.addMut(coord, s.nodeIdx, s.ch);
                    }
                }
            }
            w.flushEdges();
        }
        genomeLen = offset;
        if (tsDebug)
            std::cerr << "[debug] " << w.prefix << " blocksSeen=" << dbgSeen
                      << " blocksInReference=" << dbgExists
                      << " refLen=" << genomeLen << std::endl;
    };

    // Node table text (shared by every output tree sequence).
    auto writeNodeTable = [&](std::ostream& os) {
        os << "id\tis_sample\ttime\tmetadata\n";
        for (size_t i = 0; i < idToNode.size(); i++) {
            const panmanUtils::Node* n = tree.allNodes.at(idToNode[i]);
            const int isSample = n->children.empty() ? 1 : 0;
            os << i << '\t' << isSample << '\t' << nodeTime[idToNode[i]]
               << '\t' << idToNode[i] << '\n';
        }
    };

    auto flushWriter = [&](TableWriter& w, int64_t genomeLen) {
        {
            std::ofstream nodesOut(w.prefix + ".nodes.txt");
            writeNodeTable(nodesOut);
        }
        {
            std::ofstream edgesOut(w.prefix + ".edges.txt");
            edgesOut << "left\tright\tparent\tchild\n";
            for (const auto& e : w.edges) {
                edgesOut << std::get<0>(e) << '\t' << std::get<1>(e) << '\t'
                         << std::get<2>(e) << '\t' << std::get<3>(e) << '\n';
            }
        }
        if (emitMutations) {
            std::unordered_map<int64_t, int> siteId;
            siteId.reserve(w.sites.size());
            {
                std::ofstream sitesOut(w.prefix + ".sites.txt");
                sitesOut << "position\tancestral_state\n";
                int idx = 0;
                for (const auto& kv : w.sites) {
                    siteId[kv.first] = idx++;
                    sitesOut << kv.first << '\t' << kv.second << '\n';
                }
            }
            std::ofstream mutsOut(w.prefix + ".mutations.txt");
            mutsOut << "site\tnode\tderived_state\tparent\n";
            for (const auto& m : w.muts) {
                mutsOut << siteId[std::get<0>(m)] << '\t' << std::get<1>(m) << '\t'
                        << std::get<2>(m) << "\t-1\n";
            }
        }
        {
            std::ofstream metaOut(w.prefix + ".sequence_length.txt");
            metaOut << genomeLen << '\n';
        }
        if (emitMutations) {
            std::cout << "  wrote " << w.prefix << ".{nodes,edges,sites,mutations}.txt"
                      << "  (nodes=" << idToNode.size() << ", edges=" << w.edges.size()
                      << ", sites=" << w.sites.size() << ", mutations=" << w.muts.size()
                      << ", sequence_length=" << genomeLen << ")" << std::endl;
        } else {
            std::cout << "  wrote " << w.prefix << ".{nodes,edges}.txt"
                      << "  (nodes=" << idToNode.size() << ", edges=" << w.edges.size()
                      << ", sequence_length=" << genomeLen
                      << "; sites/mutations skipped)" << std::endl;
        }
    };

    if (perChromosome) {
        for (const auto& cs : chromBlocks) {
            TableWriter w;
            w.prefix = "./info/" + outPrefix + "." + cs.first;
            std::vector<std::pair<std::string, std::vector<int>>> single{ cs };
            int64_t genomeLen = 0;
            processChromSet(w, single, genomeLen);
            flushWriter(w, genomeLen);
        }
    } else {
        TableWriter w;
        w.prefix = "./info/" + outPrefix;
        int64_t genomeLen = 0;
        processChromSet(w, chromBlocks, genomeLen);
        flushWriter(w, genomeLen);
    }

    std::cout << "Tree-sequence tables written. Assemble .trees with:\n"
              << "  python3 scripts/panman_to_trees.py ./info/" << outPrefix
              << (perChromosome ? ".<chr>" : "") << std::endl;
}

}  // namespace panmanUtils
