// ============================================================================
//  panmanUtils::Network : DAG built from a sequence of SPR-related trees.
//
//  The first tree in the input file provides the base topology (tree 1),
//  represented with ordinary Node::parent / Node::children pointers. Each
//  subsequent tree is related to the previous one by a single SPR move, and
//  that move contributes exactly one directed "network edge" (stored
//  separately in `networkEdges`). Clade identity across trees is determined
//  by leaf-set equality.
// ============================================================================

#include <algorithm>
#include <functional>
#include <random>
#include <unordered_set>
#include <queue>

std::string panmanUtils::Network::newInternalId() {
    return "net_" + std::to_string(++m_nextInternalId);
}

panmanUtils::Network::Network(std::istream& fin) {
    // Read the file into (header, newick) pairs. A Newick line is required to
    // end with ';' (and contain parentheses); any other non-empty line that
    // precedes it is retained as that tree's header.
    std::vector< std::pair< std::string, std::string > > steps;
    std::string lastHeader;
    std::string line;
    while (std::getline(fin, line)) {
        std::string l = panmanUtils::stripString(line);
        if (l.empty()) continue;
        bool looksLikeNewick = (l.back() == ';') &&
                               (l.find('(') != std::string::npos);
        if (looksLikeNewick) {
            steps.push_back({ lastHeader, l });
            lastHeader.clear();
        } else {
            lastHeader = l;
        }
    }

    if (steps.empty()) {
        std::cerr << "Warning: Network input contained no Newick trees." << std::endl;
        return;
    }

    PNode* prevTree = nullptr;
    for (size_t k = 0; k < steps.size(); k++) {
        const std::string& header = steps[k].first;
        const std::string& newick = steps[k].second;

        PNode* curr = parseNewick(newick);
        if (!curr) {
            std::cerr << "Warning: failed to parse Newick at step " << k << std::endl;
            continue;
        }
        computeLeafSets(curr);

        if (k == 0) {
            buildTree1(curr);
        } else {
            // Parse MCC from header.
            std::set< std::string > lsM = parseMCCFromHeader(header);
            if (lsM.empty()) {
                std::cerr << "Warning: no MCC found in header for step " << k
                          << ": \"" << header << "\"" << std::endl;
            } else {
                PNode* rootM = findByLeafSet(curr, lsM);
                if (!rootM) {
                    std::cerr << "Warning: MCC {...} is not a clade in tree at step "
                              << k << ": \"" << header << "\"" << std::endl;
                } else if (!rootM->parent) {
                    std::cerr << "Warning: moved clade is the root itself at step "
                              << k << std::endl;
                } else {
                    // Find the (first) sibling of rootM in the current tree k.
                    PNode* sibling = nullptr;
                    for (PNode* c : rootM->parent->children) {
                        if (c != rootM) { sibling = c; break; }
                    }
                    if (!sibling) {
                        std::cerr << "Warning: moved clade has no sibling at step "
                                  << k << std::endl;
                    } else if (!prevTree) {
                        std::cerr << "Warning: no previous tree at step " << k << std::endl;
                    } else {
                        // Find sibling's clade in the PREVIOUS tree, then take its parent.
                        PNode* sPrev = findByLeafSet(prevTree, sibling->leafSet);
                        if (!sPrev) {
                            std::cerr << "Warning: sibling clade not present in previous "
                                      << "tree at step " << k
                                      << " (step-over-step topology mismatch)." << std::endl;
                        } else if (!sPrev->parent) {
                            std::cerr << "Warning: sibling's parent in previous tree is "
                                      << "missing at step " << k << std::endl;
                        } else {
                            const std::set< std::string >& lsPs = sPrev->parent->leafSet;
                            Node* pNode = getOrCreateNode(lsPs);
                            Node* mNode = getOrCreateNode(lsM);
                            addNetworkEdgeSafe(pNode, mNode);
                        }
                    }
                }
            }
        }

        if (prevTree) freePNode(prevTree);
        prevTree = curr;
    }
    if (prevTree) freePNode(prevTree);
}

panmanUtils::Network::~Network() {
    for (auto& kv : allNodes) {
        delete kv.second;
    }
    allNodes.clear();
}

size_t panmanUtils::Network::numLeaves() const {
    size_t n = 0;
    for (const auto& kv : allNodes) {
        if (kv.second->children.empty() && tree1Depth.count(kv.first)) n++;
    }
    return n;
}

// ----- Newick parsing -----

panmanUtils::Network::PNode* panmanUtils::Network::parseNewick(const std::string& newick) {
    std::string s = newick;
    // Strip trailing ';' and whitespace.
    while (!s.empty() && (s.back() == ';' ||
                          s.back() == ' ' || s.back() == '\t' ||
                          s.back() == '\n' || s.back() == '\r')) {
        s.pop_back();
    }
    // Strip leading whitespace.
    size_t start = 0;
    while (start < s.size() && (s[start] == ' ' || s[start] == '\t')) start++;
    s = s.substr(start);
    if (s.empty()) return nullptr;
    size_t pos = 0;
    PNode* root = parseNewickHelper(s, pos);
    return root;
}

panmanUtils::Network::PNode* panmanUtils::Network::parseNewickHelper(
    const std::string& s, size_t& pos) {
    PNode* node = new PNode();
    if (pos < s.size() && s[pos] == '(') {
        pos++; // consume '('
        while (true) {
            PNode* child = parseNewickHelper(s, pos);
            child->parent = node;
            node->children.push_back(child);
            if (pos < s.size() && s[pos] == ',') {
                pos++;
                continue;
            } else if (pos < s.size() && s[pos] == ')') {
                pos++;
                break;
            } else {
                break;
            }
        }
    }
    // Optional label.
    std::string label;
    while (pos < s.size() && s[pos] != ',' && s[pos] != ')' && s[pos] != ':') {
        label += s[pos++];
    }
    // Trim surrounding whitespace.
    while (!label.empty() && (label.front() == ' ' || label.front() == '\t'))
        label.erase(label.begin());
    while (!label.empty() && (label.back() == ' ' || label.back() == '\t'))
        label.pop_back();
    node->label = label;
    // Optional branch length (discarded).
    if (pos < s.size() && s[pos] == ':') {
        pos++;
        while (pos < s.size() && s[pos] != ',' && s[pos] != ')') pos++;
    }
    return node;
}

void panmanUtils::Network::computeLeafSets(PNode* node) {
    node->leafSet.clear();
    if (node->children.empty()) {
        if (!node->label.empty()) {
            node->leafSet.insert(node->label);
        }
    } else {
        for (PNode* c : node->children) {
            computeLeafSets(c);
            for (const std::string& l : c->leafSet) {
                node->leafSet.insert(l);
            }
        }
    }
}

void panmanUtils::Network::freePNode(PNode* node) {
    if (!node) return;
    for (PNode* c : node->children) freePNode(c);
    delete node;
}

panmanUtils::Network::PNode* panmanUtils::Network::findByLeafSet(
    PNode* root, const std::set< std::string >& ls) {
    if (!root) return nullptr;
    if (root->leafSet == ls) return root;
    for (PNode* c : root->children) {
        PNode* r = findByLeafSet(c, ls);
        if (r) return r;
    }
    return nullptr;
}

// ----- Tree-1 construction -----

void panmanUtils::Network::buildTree1(PNode* pRoot) {
    std::function< Node*(PNode*, Node*) > dfs = [&](PNode* p, Node* parNode) -> Node* {
        std::string id;
        if (p->children.empty()) {
            id = p->label.empty() ? newInternalId() : p->label;
        } else {
            id = newInternalId();
        }
        Node* n;
        if (parNode == nullptr) {
            n = new Node(id, 0.0f);
        } else {
            n = new Node(id, parNode, 0.0f);
        }
        allNodes[id] = n;
        leafsetToNode[p->leafSet] = n;
        tree1Depth[id] = n->level;
        for (PNode* c : p->children) {
            dfs(c, n);
        }
        return n;
    };
    root = dfs(pRoot, nullptr);
}

panmanUtils::Node* panmanUtils::Network::getOrCreateNode(
    const std::set< std::string >& ls) {
    auto it = leafsetToNode.find(ls);
    if (it != leafsetToNode.end()) return it->second;
    std::string id;
    if (ls.size() == 1) {
        // A singleton clade should correspond to a leaf that already exists if
        // it was part of tree 1. If we get here, the leaf did not appear in
        // tree 1 -- create a fresh leaf node for it.
        id = *ls.begin();
    } else {
        id = newInternalId();
    }
    Node* n = new Node(id, 0.0f);
    allNodes[id] = n;
    leafsetToNode[ls] = n;
    return n;
}

std::set< std::string > panmanUtils::Network::parseMCCFromHeader(
    const std::string& header) {
    std::set< std::string > result;
    size_t mccPos = header.find("MCC");
    if (mccPos == std::string::npos) return result;
    size_t lb = header.find('{', mccPos);
    size_t rb = (lb == std::string::npos) ? std::string::npos : header.find('}', lb);
    if (lb == std::string::npos || rb == std::string::npos) return result;
    std::string content = header.substr(lb + 1, rb - lb - 1);
    std::string cur;
    auto push = [&](std::string tok) {
        tok = panmanUtils::stripString(tok);
        if (!tok.empty()) result.insert(tok);
    };
    for (char c : content) {
        if (c == ',') {
            push(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    push(cur);
    return result;
}

// ----- Network edge addition & conflict resolution -----

bool panmanUtils::Network::isAncestor(Node* anc, Node* desc) const {
    if (!anc || !desc || anc == desc) return anc == desc;
    std::queue< Node* > q;
    std::unordered_set< Node* > visited;
    q.push(desc);
    visited.insert(desc);
    while (!q.empty()) {
        Node* n = q.front(); q.pop();
        if (n == anc) return true;
        // Follow tree-parent pointer.
        if (n->parent && !visited.count(n->parent)) {
            visited.insert(n->parent);
            q.push(n->parent);
        }
        // Follow network edges where `n` is the child.
        for (const NetworkEdge& e : networkEdges) {
            if (e.child == n && !visited.count(e.parent)) {
                visited.insert(e.parent);
                q.push(e.parent);
            }
        }
    }
    return false;
}

void panmanUtils::Network::addNetworkEdgeSafe(Node* par, Node* ch) {
    if (!par || !ch || par == ch) return;
    // If par is already the tree parent of ch, nothing to do.
    if (ch->parent == par) return;
    // Deduplicate network edges.
    for (const NetworkEdge& e : networkEdges) {
        if (e.parent == par && e.child == ch) return;
    }
    // Cycle check: is `ch` already an ancestor of `par`?
    if (isAncestor(ch, par)) {
        auto getDepth = [&](Node* n) -> size_t {
            auto it = tree1Depth.find(n->identifier);
            return (it == tree1Depth.end()) ? std::numeric_limits< size_t >::max()
                                            : it->second;
        };
        size_t pd = getDepth(par);
        size_t cd = getDepth(ch);
        // Rule: the deeper (tree-1) node stays as the child. The incoming edge
        // would make `ch` a child of `par`, so accepting it requires that `ch`
        // be deeper (>=) than `par` in tree 1. Otherwise we reject the
        // incoming edge.
        if (cd >= pd) {
            std::cerr << "Warning: SPR edge " << par->identifier << " -> "
                      << ch->identifier
                      << " rejected (would create a cycle; existing tree-1 "
                         "ancestry is preferred)." << std::endl;
            return;
        } else {
            // `par` is deeper than `ch` in tree 1, so `par` should be the
            // child. The requested edge (par -> ch) contradicts that, so we
            // still reject it; future extensions may remove the conflicting
            // ancestry instead.
            std::cerr << "Warning: SPR edge " << par->identifier << " -> "
                      << ch->identifier
                      << " rejected (par is deeper than ch in tree 1; conflict)."
                      << std::endl;
            return;
        }
    }
    networkEdges.push_back({ par, ch });
}

// ----- Summary / printing -----

void panmanUtils::Network::printSummary(std::ostream& out) const {
    size_t nLeaves = 0, nInternal = 0, nExtra = 0;
    size_t treeEdges = 0;
    for (const auto& kv : allNodes) {
        Node* n = kv.second;
        bool inTree1 = tree1Depth.count(n->identifier) > 0;
        if (!inTree1) {
            nExtra++;
        } else if (n->children.empty()) {
            nLeaves++;
        } else {
            nInternal++;
        }
        if (n->parent) treeEdges++;
    }
    out << "Network summary" << std::endl;
    out << "  Tree-1 leaves              : " << nLeaves << std::endl;
    out << "  Tree-1 internal nodes      : " << nInternal << std::endl;
    out << "  Extra network-only nodes   : " << nExtra << std::endl;
    out << "  Tree-1 parent->child edges : " << treeEdges << std::endl;
    out << "  SPR network edges          : " << networkEdges.size() << std::endl;
}

void panmanUtils::Network::printNewick(std::ostream& out, Node* node) const {
    if (!node) return;
    if (!node->children.empty()) {
        out << '(';
        for (size_t i = 0; i < node->children.size(); i++) {
            if (i) out << ',';
            printNewick(out, node->children[i]);
        }
        out << ')';
    }
    out << node->identifier;
}

void panmanUtils::Network::printNetwork(std::ostream& out) const {
    out << "# base tree (tree 1) in Newick format" << std::endl;
    if (root) {
        printNewick(out, root);
        out << ';' << std::endl;
    }
    out << std::endl;
    out << "# network edges (parent -> child), added from SPR moves" << std::endl;
    for (const NetworkEdge& e : networkEdges) {
        out << e.parent->identifier << " -> " << e.child->identifier << std::endl;
    }
}

// ============================================================================
//  buildNetworkFromTreeknit : construct a PanMAN (TreeGroup) from TreeKnit output.
//
//  Layout (see panman.hpp declaration): one PanMAT per window, each window a
//  single block (uncovered reference gaps merged into the left neighbour block),
//  adjacent windows linked by one SPR-derived recombination complex mutation.
// ============================================================================

#include <filesystem>
#include <algorithm>
#include <fstream>
#include <cctype>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/lzma.hpp>

namespace {

namespace fs = std::filesystem;

struct TkWindow {
    std::string name;        // chrXX_A_B
    std::string chrom;       // chrXX
    long long start = 0;     // reference start (A)
    long long end = 0;       // reference end (B)
    std::string newickPath;  // path to <name>_resolved.nwk or <name>.nwk
    int chrIdx = 0;          // chromosome index (for block annotation)
};

struct TkMove {
    size_t fromIdx = 0;                  // earlier (parent) window
    size_t toIdx = 0;                    // later (child) window
    std::vector< std::string > movedClade; // smaller MCC leaf set
};

bool tkEndsWith(const std::string& s, const std::string& suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Parse "chrXX_A_B" -> chrom="chrXX", a=A, b=B (last two underscore tokens are numeric).
bool parseWindowName(const std::string& name, std::string& chrom,
                     long long& a, long long& b) {
    std::vector< std::string > toks;
    size_t prev = 0, pos;
    while ((pos = name.find('_', prev)) != std::string::npos) {
        toks.push_back(name.substr(prev, pos - prev));
        prev = pos + 1;
    }
    toks.push_back(name.substr(prev));
    if (toks.size() < 3) return false;
    try {
        b = std::stoll(toks[toks.size() - 1]);
        a = std::stoll(toks[toks.size() - 2]);
    } catch (...) {
        return false;
    }
    chrom.clear();
    for (size_t i = 0; i + 2 < toks.size(); i++) {
        if (i) chrom += "_";
        chrom += toks[i];
    }
    return true;
}

// Read the (ungapped, upper-case) reference sequence for `chrom` from a FASTA that
// may be plain, gzip- or xz/lzma-compressed. Returns "" if not found.
std::string loadReferenceChromosome(const std::string& refFile, const std::string& chrom) {
    std::ifstream raw(refFile, std::ios::binary);
    if (!raw.is_open()) {
        std::cerr << "Warning: could not open reference file " << refFile
                  << "; skipping gap merge." << std::endl;
        return "";
    }
    boost::iostreams::filtering_streambuf< boost::iostreams::input > inbuf;
    if (tkEndsWith(refFile, ".gz")) {
        inbuf.push(boost::iostreams::gzip_decompressor());
    } else if (tkEndsWith(refFile, ".xz") || tkEndsWith(refFile, ".lzma")) {
        inbuf.push(boost::iostreams::lzma_decompressor());
    }
    inbuf.push(raw);
    std::istream in(&inbuf);

    auto matches = [&](const std::string& id) {
        if (id == chrom) return true;
        if (("chr" + id) == chrom) return true;
        if (id == ("chr" + chrom)) return true;
        if (chrom.rfind("chr", 0) == 0 && id == chrom.substr(3)) return true;
        return false;
    };

    std::string line, seq;
    bool inWanted = false;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') {
            if (inWanted) break; // finished the wanted record
            std::string id = line.substr(1);
            size_t sp = id.find_first_of(" \t");
            if (sp != std::string::npos) id = id.substr(0, sp);
            inWanted = matches(id);
        } else if (inWanted) {
            seq += line;
        }
    }
    for (auto& c : seq) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    return seq;
}

// Copy a FASTA MSA, appending `suffix` (invariant reference columns) to every record.
bool writeAugmentedMSA(const std::string& inPath, const std::string& suffix,
                       const std::string& outPath) {
    std::ifstream in(inPath);
    if (!in.is_open()) return false;
    std::ofstream out(outPath);
    if (!out.is_open()) return false;

    std::string line, curId, curSeq;
    bool have = false;
    auto flush = [&]() {
        if (have) out << ">" << curId << "\n" << curSeq << suffix << "\n";
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') {
            flush();
            curId = line.substr(1);
            curSeq.clear();
            have = true;
        } else {
            curSeq += line;
        }
    }
    flush();
    return true;
}

// Most recent common ancestor of the given leaf set within tree T (by Node::level).
panmanUtils::Node* findMRCA(panmanUtils::Tree& T,
                            const std::vector< std::string >& leaves) {
    panmanUtils::Node* mrca = nullptr;
    for (const std::string& ln : leaves) {
        auto it = T.allNodes.find(ln);
        if (it == T.allNodes.end()) continue;
        panmanUtils::Node* n = it->second;
        if (!mrca) {
            mrca = n;
            continue;
        }
        panmanUtils::Node* a = mrca;
        panmanUtils::Node* b = n;
        while (a->level > b->level) a = a->parent;
        while (b->level > a->level) b = b->parent;
        while (a != b && a && b) {
            a = a->parent;
            b = b->parent;
        }
        if (a) mrca = a;
    }
    return mrca;
}

} // namespace

panmanUtils::TreeGroup* panmanUtils::buildNetworkFromTreeknit(
    const std::string& treeknitDir, const std::string& alignmentDir,
    const std::string& refFile) {

    // ----- 1. Enumerate window directories -----
    std::vector< fs::path > windowDirs;
    std::error_code ec;
    if (!fs::exists(treeknitDir, ec)) {
        std::cerr << "Error: TreeKnit directory does not exist: " << treeknitDir << std::endl;
        return nullptr;
    }
    for (const auto& e : fs::directory_iterator(treeknitDir, ec)) {
        if (e.is_directory(ec) &&
            e.path().filename().string().find("_window_") != std::string::npos) {
            windowDirs.push_back(e.path());
        }
    }
    if (windowDirs.empty()) {
        std::cerr << "Error: no chrXX_window_N_N+1 subdirectories found in "
                  << treeknitDir << std::endl;
        return nullptr;
    }
    std::sort(windowDirs.begin(), windowDirs.end());

    // ----- 2. Parse MCCs.json in each window directory -----
    struct RawMove {
        std::string treeA, treeB;
        std::vector< std::string > movedClade;
    };
    std::map< std::string, TkWindow > windowMap;
    std::vector< RawMove > rawMoves;

    for (const auto& d : windowDirs) {
        fs::path mccPath = d / "MCCs.json";
        std::ifstream jf(mccPath);
        if (!jf.is_open()) {
            std::cerr << "Warning: missing " << mccPath << std::endl;
            continue;
        }
        Json::Value jroot;
        try {
            jf >> jroot;
        } catch (const std::exception& e) {
            std::cerr << "Warning: failed to parse " << mccPath << ": " << e.what()
                      << std::endl;
            continue;
        }
        const Json::Value& dict = jroot["MCC_dict"];
        for (const auto& key : dict.getMemberNames()) {
            const Json::Value& entry = dict[key];
            if (!entry.isMember("trees") || entry["trees"].size() < 2) continue;
            std::string ta = entry["trees"][0].asString();
            std::string tb = entry["trees"][1].asString();

            const Json::Value& mccs = entry["mccs"];
            if (mccs.size() < 2) continue;
            std::vector< std::string > m0, m1;
            for (const auto& x : mccs[0]) m0.push_back(x.asString());
            for (const auto& x : mccs[1]) m1.push_back(x.asString());
            std::vector< std::string > moved = (m0.size() <= m1.size()) ? m0 : m1;

            for (const std::string& tn : { ta, tb }) {
                std::string chrom;
                long long a = 0, b = 0;
                if (!parseWindowName(tn, chrom, a, b)) continue;
                auto it = windowMap.find(tn);
                if (it == windowMap.end()) {
                    TkWindow w;
                    w.name = tn;
                    w.chrom = chrom;
                    w.start = a;
                    w.end = b;
                    windowMap[tn] = w;
                    it = windowMap.find(tn);
                }
                if (it->second.newickPath.empty()) {
                    fs::path nwk = d / (tn + "_resolved.nwk");
                    if (fs::exists(nwk, ec)) it->second.newickPath = nwk.string();
                }
            }
            rawMoves.push_back({ ta, tb, moved });
        }
    }

    if (windowMap.empty()) {
        std::cerr << "Error: no valid windows parsed from MCCs.json files." << std::endl;
        return nullptr;
    }

    // ----- 3. Order windows by reference start coordinate -----
    std::vector< TkWindow > windows;
    windows.reserve(windowMap.size());
    for (auto& kv : windowMap) windows.push_back(kv.second);
    std::sort(windows.begin(), windows.end(),
              [](const TkWindow& x, const TkWindow& y) { return x.start < y.start; });

    std::map< std::string, size_t > nameToIdx;
    for (size_t i = 0; i < windows.size(); i++) nameToIdx[windows[i].name] = i;

    // Fill in any missing newick paths by searching all window directories.
    for (auto& w : windows) {
        if (!w.newickPath.empty()) continue;
        for (const auto& d : windowDirs) {
            fs::path nwk = d / (w.name + "_resolved.nwk");
            if (fs::exists(nwk, ec)) {
                w.newickPath = nwk.string();
                break;
            }
        }
    }

    // ----- 4. Convert raw moves to ordered (fromIdx<toIdx) moves -----
    std::vector< TkMove > moves;
    for (const auto& rm : rawMoves) {
        if (!nameToIdx.count(rm.treeA) || !nameToIdx.count(rm.treeB)) continue;
        size_t ia = nameToIdx[rm.treeA];
        size_t ib = nameToIdx[rm.treeB];
        TkMove mv;
        mv.fromIdx = std::min(ia, ib);
        mv.toIdx = std::max(ia, ib);
        mv.movedClade = rm.movedClade;
        moves.push_back(mv);
    }
    std::sort(moves.begin(), moves.end(),
              [](const TkMove& x, const TkMove& y) { return x.toIdx < y.toIdx; });

    // ----- 5. Load reference chromosome (single chromosome assumed across windows) -----
    std::string refSeq = refFile.empty() ? std::string()
                                         : loadReferenceChromosome(refFile, windows[0].chrom);

    // ----- 6. Build one PanMAT per window (FILE_TYPE::MSA, reference GRCh38) -----
    fs::path tmpDir = fs::path(".") / "tmp_treeknit_msa";
    fs::create_directories(tmpDir, ec);

    std::vector< panmanUtils::Tree* > treePtrs(windows.size(), nullptr);
    for (size_t i = 0; i < windows.size(); i++) {
        const TkWindow& w = windows[i];
        std::string msaPath = (fs::path(alignmentDir) / (w.name + ".fa")).string();
        if (!fs::exists(msaPath, ec)) {
            std::cerr << "Error: missing alignment file " << msaPath << std::endl;
            return nullptr;
        }
        if (w.newickPath.empty()) {
            std::cerr << "Error: missing resolved newick for window " << w.name << std::endl;
            return nullptr;
        }

        // Gap merge: append uncovered reference region [end, nextStart) to the LEFT
        // (current) window's block as invariant columns.
        std::string suffix;
        if (!refSeq.empty() && i + 1 < windows.size()) {
            long long gapBeg = w.end;                  // 0-based index of first uncovered base
            long long gapLen = windows[i + 1].start - w.end - 1;
            if (gapLen > 0 && gapBeg >= 0 &&
                gapBeg + gapLen <= static_cast<long long>(refSeq.size())) {
                suffix = refSeq.substr(static_cast<size_t>(gapBeg),
                                       static_cast<size_t>(gapLen));
            }
        }

        std::string useMsaPath = msaPath;
        std::string tmpMsa;
        if (!suffix.empty()) {
            tmpMsa = (tmpDir / (w.name + ".aug.fa")).string();
            if (writeAugmentedMSA(msaPath, suffix, tmpMsa)) {
                useMsaPath = tmpMsa;
            } else {
                std::cerr << "Warning: gap merge failed for " << w.name
                          << "; using original MSA." << std::endl;
                tmpMsa.clear();
            }
        }

        std::ifstream msaIfs(useMsaPath);
        std::ifstream nwkIfs(w.newickPath);
        std::cout << "Building PanMAT for window " << w.name << " (" << (i + 1) << "/"
                  << windows.size() << ")" << std::endl;
        panmanUtils::Tree* T = new panmanUtils::Tree(
            msaIfs, nwkIfs, panmanUtils::FILE_TYPE::MSA, std::string("GRCh38"));
        treePtrs[i] = T;
        msaIfs.close();
        nwkIfs.close();
        if (!tmpMsa.empty()) fs::remove(tmpMsa, ec);
    }

    // ----- 7. Assemble TreeGroup -----
    panmanUtils::TreeGroup* TG = new panmanUtils::TreeGroup(treePtrs);
    if (!TG->trees.empty() && TG->trees[0].root) {
        // Base tree (window 0) is the head/pseudo-root of the network.
        TG->trees[0].root->isComMutHead = true;
        TG->trees[0].root->treeIndex = 0;
    }

    // ----- 8. One recombination complex mutation per SPR move -----
    for (const auto& mv : moves) {
        if (mv.toIdx >= TG->trees.size() || mv.fromIdx >= TG->trees.size()) continue;
        panmanUtils::Tree& childTree = TG->trees[mv.toIdx];
        panmanUtils::Tree& parentTree = TG->trees[mv.fromIdx];

        panmanUtils::Node* childMRCA = findMRCA(childTree, mv.movedClade);
        panmanUtils::Node* parentMRCA = findMRCA(parentTree, mv.movedClade);
        if (!childMRCA || !parentMRCA) {
            std::cerr << "Warning: could not locate moved clade for move "
                      << windows[mv.fromIdx].name << " -> " << windows[mv.toIdx].name
                      << "; skipping." << std::endl;
            continue;
        }

        // Breakpoints: the parent's window block (MSA coordinates) contributes the
        // recombinant segment. Single block per window -> primaryBlockId 0.
        long long parentBlockLen =
            parentTree.blocks.empty() ? 0LL
                                      : static_cast<long long>(parentTree.blocks[0].blockLength);
        int lastNuc = static_cast<int>(parentBlockLen > 0 ? parentBlockLen - 1 : 0);
        std::tuple< int, int, int, int > startCoords(0, -1, 0, -1);
        std::tuple< int, int, int, int > endCoords(0, -1, lastNuc, -1);

        std::vector< panmanUtils::ParentContribution > parents;
        parents.emplace_back(mv.fromIdx, parentMRCA->identifier, startCoords, endCoords);
        TG->complexMutations.emplace_back('R', mv.toIdx, childMRCA->identifier,
                                          std::move(parents));

        // Pseudo-root bookkeeping: this child tree has a parent-tree connection.
        if (childTree.root) {
            childTree.root->isComMutHead = true;
            childTree.root->treeIndex = static_cast<int>(mv.toIdx);
        }
    }

    fs::remove_all(tmpDir, ec);

    std::cout << "Built network from TreeKnit: " << TG->trees.size() << " trees, "
              << TG->complexMutations.size() << " complex mutations." << std::endl;
    return TG;
}

// ============================================================================
//  buildNetworkFromSprDirs : single-tree + network-edges architecture
// ============================================================================

namespace {

// ----- debug helpers -----
bool sprDebugEnabled() {
    static int val = -1;
    if (val < 0) val = (std::getenv("PANMAN_SPR_DEBUG") != nullptr) ? 1 : 0;
    return val == 1;
}
void sprDebugLog(const std::string& msg) {
    std::cerr << "[SPR-DBG] " << msg << std::endl;
}
void sprDebugClade(const char* label, const std::vector< std::string >& clade) {
    std::ostringstream os;
    os << "  " << label << " = {";
    for (size_t i = 0; i < clade.size(); i++) {
        if (i) os << ", ";
        os << clade[i];
    }
    os << "}";
    sprDebugLog(os.str());
}
void sprDebugTreeStats(const char* label, const panmanUtils::Tree& T) {
    std::ostringstream os;
    os << label << ": nodes=" << T.allNodes.size()
       << " blocks=" << T.blocks.size()
       << " gaps=" << T.gaps.size();
    sprDebugLog(os.str());
}

// ----- MSA path resolution -----
std::string resolveMsaPath(const std::string& msaRoot, const TkWindow& w) {
    fs::path nested = fs::path(msaRoot) / w.chrom / "full_msa" / (w.name + ".fa");
    if (fs::exists(nested)) return nested.string();
    fs::path flat = fs::path(msaRoot) / (w.chrom + "_full_msa") / (w.name + ".fa");
    if (fs::exists(flat)) return flat.string();
    fs::path nestedFasta = fs::path(msaRoot) / w.chrom / "full_msa" / (w.name + ".fasta");
    if (fs::exists(nestedFasta)) return nestedFasta.string();
    return "";
}

// ----- Read a FASTA MSA into a map<seqName, sequence> -----
size_t readMsaFile(const std::string& path,
                   std::map< std::string, std::string >& seqs) {
    seqs.clear();
    std::ifstream in(path);
    if (!in.is_open()) return 0;
    std::string line, curId, curSeq;
    auto flush = [&]() {
        if (!curId.empty() && !curSeq.empty()) {
            seqs[curId] = curSeq;
        }
    };
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == '>') {
            flush();
            curId = line.substr(1);
            size_t sp = curId.find_first_of(" \t");
            if (sp != std::string::npos) curId = curId.substr(0, sp);
            curSeq.clear();
        } else {
            curSeq += line;
        }
    }
    flush();
    size_t alignLen = 0;
    for (const auto& kv : seqs) {
        if (kv.second.size() > alignLen) alignLen = kv.second.size();
    }
    return alignLen;
}

// ----- SPR inference from consecutive Newick trees -----
struct SprMoveInfo {
    std::vector< std::string > movedClade;
    std::vector< std::string > acceptorClade;
    std::vector< std::string > donorClade;
};

struct SprStep {
    size_t windowIdx;
    SprMoveInfo spr;
};

struct SimpleNode {
    std::string label;
    SimpleNode* parent = nullptr;
    std::vector< SimpleNode* > children;
    std::set< std::string > leafSet;
};

SimpleNode* parseSimpleNewick(const std::string& newick) {
    std::string s = newick;
    while (!s.empty() && (s.back() == ';' || s.back() == ' ' ||
                          s.back() == '\t' || s.back() == '\n' ||
                          s.back() == '\r')) s.pop_back();
    size_t start = 0;
    while (start < s.size() && (s[start] == ' ' || s[start] == '\t')) start++;
    s = s.substr(start);
    if (s.empty()) return nullptr;

    std::function< SimpleNode*(const std::string&, size_t&) > parse =
        [&](const std::string& str, size_t& pos) -> SimpleNode* {
        SimpleNode* node = new SimpleNode();
        if (pos < str.size() && str[pos] == '(') {
            pos++;
            while (true) {
                SimpleNode* child = parse(str, pos);
                child->parent = node;
                node->children.push_back(child);
                if (pos < str.size() && str[pos] == ',') { pos++; continue; }
                if (pos < str.size() && str[pos] == ')') { pos++; break; }
                break;
            }
        }
        std::string label;
        while (pos < str.size() && str[pos] != ',' && str[pos] != ')' &&
               str[pos] != ':') {
            label += str[pos++];
        }
        while (!label.empty() && (label.front() == ' ' || label.front() == '\t'))
            label.erase(label.begin());
        while (!label.empty() && (label.back() == ' ' || label.back() == '\t'))
            label.pop_back();
        node->label = label;
        if (pos < str.size() && str[pos] == ':') {
            pos++;
            while (pos < str.size() && str[pos] != ',' && str[pos] != ')') pos++;
        }
        return node;
    };
    size_t pos = 0;
    return parse(s, pos);
}

void computeSimpleLeafSets(SimpleNode* node) {
    node->leafSet.clear();
    if (node->children.empty()) {
        if (!node->label.empty()) node->leafSet.insert(node->label);
    } else {
        for (auto* c : node->children) {
            computeSimpleLeafSets(c);
            node->leafSet.insert(c->leafSet.begin(), c->leafSet.end());
        }
    }
}

void freeSimpleNode(SimpleNode* node) {
    if (!node) return;
    for (auto* c : node->children) freeSimpleNode(c);
    delete node;
}

SimpleNode* findSimpleByLeafSet(SimpleNode* root,
                                const std::set< std::string >& ls) {
    if (!root) return nullptr;
    if (root->leafSet == ls) return root;
    for (auto* c : root->children) {
        auto* r = findSimpleByLeafSet(c, ls);
        if (r) return r;
    }
    return nullptr;
}

bool resolveSprMove(SimpleNode* pRoot, SimpleNode* movedNode, SprMoveInfo& info) {
    if (!movedNode || !movedNode->parent) return false;

    SimpleNode* mPrev = findSimpleByLeafSet(pRoot, movedNode->leafSet);
    if (!mPrev || !mPrev->parent) return false;

    for (SimpleNode* sibling : movedNode->parent->children) {
        if (sibling == movedNode) continue;
        SimpleNode* sPrev = findSimpleByLeafSet(pRoot, sibling->leafSet);
        if (!sPrev || !sPrev->parent) continue;

        info.movedClade.assign(movedNode->leafSet.begin(), movedNode->leafSet.end());
        info.acceptorClade.assign(sPrev->parent->leafSet.begin(),
                                  sPrev->parent->leafSet.end());
        info.donorClade.assign(mPrev->parent->leafSet.begin(),
                               mPrev->parent->leafSet.end());
        return true;
    }
    return false;
}

bool inferSprMove(const std::string& prevNewick, const std::string& currNewick,
                  SprMoveInfo& spr) {
    SimpleNode* pRoot = parseSimpleNewick(prevNewick);
    SimpleNode* cRoot = parseSimpleNewick(currNewick);
    if (!pRoot || !cRoot) {
        freeSimpleNode(pRoot);
        freeSimpleNode(cRoot);
        return false;
    }
    computeSimpleLeafSets(pRoot);
    computeSimpleLeafSets(cRoot);

    std::vector< SimpleNode* > candidates;
    std::function< void(SimpleNode*) > walk = [&](SimpleNode* n) {
        if (!n->leafSet.empty() && n->parent) {
            SimpleNode* pMatch = findSimpleByLeafSet(pRoot, n->leafSet);
            if (pMatch && pMatch->parent && n->parent &&
                pMatch->parent->leafSet != n->parent->leafSet) {
                candidates.push_back(n);
            }
        }
        for (SimpleNode* c : n->children) walk(c);
    };
    walk(cRoot);

    if (candidates.empty()) {
        freeSimpleNode(pRoot);
        freeSimpleNode(cRoot);
        return false;
    }

    std::sort(candidates.begin(), candidates.end(),
              [](SimpleNode* a, SimpleNode* b) {
                  if (a->leafSet.size() != b->leafSet.size()) {
                      return a->leafSet.size() < b->leafSet.size();
                  }
                  return a->leafSet < b->leafSet;
              });

    bool found = false;
    size_t minSize = candidates.front()->leafSet.size();

    for (SimpleNode* movedNode : candidates) {
        if (movedNode->leafSet.size() > minSize) break;
        SprMoveInfo attempt;
        if (resolveSprMove(pRoot, movedNode, attempt)) {
            spr = std::move(attempt);
            found = true;
            break;
        }
    }
    if (!found) {
        for (SimpleNode* movedNode : candidates) {
            if (movedNode->leafSet.size() == minSize) continue;
            SprMoveInfo attempt;
            if (resolveSprMove(pRoot, movedNode, attempt)) {
                spr = std::move(attempt);
                found = true;
                break;
            }
        }
    }

    freeSimpleNode(pRoot);
    freeSimpleNode(cRoot);
    return found;
}

std::vector< SprStep > collectSprSteps(const std::vector< TkWindow >& windows) {
    std::vector< SprStep > steps;
    for (size_t i = 1; i < windows.size(); i++) {
        if (windows[i].newickPath.empty() || windows[i - 1].newickPath.empty())
            continue;
        std::ifstream prevIfs(windows[i - 1].newickPath);
        std::ifstream currIfs(windows[i].newickPath);
        std::string prevNwk, currNwk;
        std::getline(prevIfs, prevNwk);
        std::getline(currIfs, currNwk);

        SprMoveInfo spr;
        if (prevNwk == currNwk) {
            if (sprDebugEnabled()) {
                sprDebugLog("identical trees for "
                            + windows[i - 1].name + " -> " + windows[i].name
                            + "; no SPR.");
            }
        } else if (inferSprMove(prevNwk, currNwk, spr)) {
            steps.push_back({ i, spr });
        } else {
            std::cerr << "Warning: could not infer SPR move for "
                      << windows[i - 1].name << " -> " << windows[i].name
                      << "; skipping." << std::endl;
        }
    }
    return steps;
}

// ----- SPR application/undo on tree -----
struct AppliedSpr {
    panmanUtils::Node* node;
    panmanUtils::Node* originalParent;
    size_t originalChildIndex;
};

AppliedSpr applySprToTree(panmanUtils::Tree& tree,
                          panmanUtils::Node* movedNode,
                          panmanUtils::Node* newParent) {
    AppliedSpr result;
    result.node = movedNode;
    result.originalParent = movedNode->parent;
    result.originalChildIndex = 0;
    if (movedNode->parent) {
        auto& siblings = movedNode->parent->children;
        for (size_t i = 0; i < siblings.size(); i++) {
            if (siblings[i] == movedNode) {
                result.originalChildIndex = i;
                siblings.erase(siblings.begin() + i);
                break;
            }
        }
    }
    movedNode->parent = newParent;
    newParent->children.push_back(movedNode);
    return result;
}

void undoSpr(const AppliedSpr& applied) {
    if (!applied.node || !applied.originalParent) return;
    panmanUtils::Node* node = applied.node;
    panmanUtils::Node* currParent = node->parent;
    if (currParent) {
        auto& siblings = currParent->children;
        for (auto it = siblings.begin(); it != siblings.end(); ++it) {
            if (*it == node) {
                siblings.erase(it);
                break;
            }
        }
    }
    node->parent = applied.originalParent;
    auto& targetChildren = applied.originalParent->children;
    size_t idx = std::min(applied.originalChildIndex, targetChildren.size());
    targetChildren.insert(targetChildren.begin() + idx, node);
}

bool isAncestorOf(panmanUtils::Node* anc, panmanUtils::Node* desc) {
    if (!anc || !desc) return false;
    panmanUtils::Node* cur = desc;
    while (cur) {
        if (cur == anc) return true;
        cur = cur->parent;
    }
    return false;
}

// ----- Run Fitch for a single block on the current tree topology -----
void runFitchForBlock(panmanUtils::Tree& tree,
                      const std::map< std::string, std::string >& leafSeqs,
                      size_t alignLen, int32_t blockId,
                      const std::string& refName) {
    panmanUtils::Alphabet alphabet = panmanUtils::getActiveAlphabet();
    std::string consensusSeq(alignLen, 'N');

    // Find reference sequence if available
    std::string refSeqStr;
    if (!refName.empty()) {
        auto it = leafSeqs.find(refName);
        if (it != leafSeqs.end()) refSeqStr = it->second;
    }
    if (!refSeqStr.empty()) {
        consensusSeq = refSeqStr;
    }

    tbb::concurrent_unordered_map< std::string,
        std::vector< std::tuple< int, int8_t, int8_t > > > nonGapMutations;
    std::unordered_map< std::string, std::mutex > nodeMutexes;
    for (const auto& kv : tree.allNodes) {
        nodeMutexes[kv.first];
    }

    tbb::parallel_for((size_t)0, alignLen, [&](size_t i) {
        std::unordered_map< std::string, int > states;
        std::unordered_map< std::string,
            std::pair< panmanUtils::NucMutationType, char > > mutations;

        for (const auto& kv : leafSeqs) {
            char c = (i < kv.second.size()) ? kv.second[i] : '-';
            if (c != '-') {
                states[kv.first] = (1 << panmanUtils::getCodeFromSymbol(c, alphabet));
            } else {
                states[kv.first] = 1;
            }
        }

        int refState = -1;
        if (!refSeqStr.empty() && i < refSeqStr.size()) {
            refState = 1 << panmanUtils::getCodeFromSymbol(refSeqStr[i], alphabet);
        }

        tree.nucFitchForwardPass(tree.root, states, refState);

        if (states[tree.root->identifier] != 1) {
            int codeLocal = 0, currentState = states[tree.root->identifier];
            while (currentState > 0) { currentState >>= 1; codeLocal++; }
            codeLocal--;
            consensusSeq[i] = panmanUtils::getSymbolFromCode(codeLocal, alphabet);
        } else {
            consensusSeq[i] = '-';
        }

        int rootState = (1 << panmanUtils::getCodeFromSymbol(consensusSeq[i], alphabet));
        tree.nucFitchBackwardPass(tree.root, states, rootState);
        tree.nucFitchAssignMutations(tree.root, states, mutations, rootState);

        for (const auto& mut : mutations) {
            nodeMutexes[mut.first].lock();
            nonGapMutations[mut.first].push_back(std::make_tuple(
                (int)i, (int8_t)mut.second.first,
                (int8_t)panmanUtils::getCodeFromSymbol(mut.second.second, alphabet)));
            nodeMutexes[mut.first].unlock();
        }
    });

    tree.blocks.emplace_back((size_t)blockId, consensusSeq, alphabet, (int64_t)alignLen);
    tree.root->blockMutation.emplace_back(
        (size_t)blockId, std::make_pair(panmanUtils::BlockMutationType::BI, false));

    tbb::parallel_for_each(nonGapMutations, [&](auto& u) {
        nodeMutexes[u.first].lock();
        std::sort(u.second.begin(), u.second.end());
        nodeMutexes[u.first].unlock();
        size_t currentStart = 0;
        for (size_t j = 1; j < u.second.size(); j++) {
            if (j - currentStart == panmanUtils::mutationPayloadCapacity(alphabet) ||
                std::get<0>(u.second[j]) != std::get<0>(u.second[j - 1]) + 1 ||
                std::get<1>(u.second[j]) != std::get<1>(u.second[j - 1])) {
                nodeMutexes[u.first].lock();
                tree.allNodes[u.first]->nucMutation.emplace_back(u.second, currentStart, j);
                nodeMutexes[u.first].unlock();
                currentStart = j;
            }
        }
        nodeMutexes[u.first].lock();
        tree.allNodes[u.first]->nucMutation.emplace_back(
            u.second, currentStart, u.second.size());
        nodeMutexes[u.first].unlock();
    });
}

// ----- Build a Tree from MSA + Newick (for window 0) -----
panmanUtils::Tree* buildWindowPanmat(
    const TkWindow& w,
    const std::vector< TkWindow >& windows,
    size_t wIdx,
    const std::string& msaRoot,
    const fs::path& tmpDir,
    std::function< const std::string&(const std::string&) > getRefSeq) {

    std::string msaPath = resolveMsaPath(msaRoot, w);
    if (msaPath.empty()) {
        std::cerr << "Error: no MSA found for " << w.name << std::endl;
        return nullptr;
    }
    if (w.newickPath.empty()) {
        std::cerr << "Error: no Newick found for " << w.name << std::endl;
        return nullptr;
    }

    std::string suffix;
    if (wIdx + 1 < windows.size() && windows[wIdx + 1].chrom == w.chrom) {
        const std::string& refSeq = getRefSeq(w.chrom);
        if (!refSeq.empty()) {
            long long gapBeg = w.end;
            long long gapLen = windows[wIdx + 1].start - w.end - 1;
            if (gapLen > 0 && gapBeg >= 0 &&
                gapBeg + gapLen <= static_cast<long long>(refSeq.size())) {
                suffix = refSeq.substr(static_cast<size_t>(gapBeg),
                                       static_cast<size_t>(gapLen));
            }
        }
    }

    std::string useMsaPath = msaPath;
    std::string tmpMsa;
    if (!suffix.empty()) {
        tmpMsa = (tmpDir / (w.name + ".aug.fa")).string();
        if (writeAugmentedMSA(msaPath, suffix, tmpMsa)) {
            useMsaPath = tmpMsa;
        } else {
            tmpMsa.clear();
        }
    }

    std::ifstream msaIfs(useMsaPath);
    std::ifstream nwkIfs(w.newickPath);
    panmanUtils::Tree* T = new panmanUtils::Tree(
        msaIfs, nwkIfs, panmanUtils::FILE_TYPE::MSA, std::string("GRCh38"));
    msaIfs.close();
    nwkIfs.close();

    std::error_code ec;
    if (!tmpMsa.empty()) fs::remove(tmpMsa, ec);
    return T;
}

// ----- Block chromosome annotation -----
void annotateBlockChromosome(panmanUtils::Tree& tree, int chrIdx,
                             const std::string& chromName,
                             int32_t blockId = 0) {
    int64_t blockIdEncoded = ((int64_t)blockId << 32);
    bool found = false;
    for (auto& chr : tree.chrList) {
        if (chr.chrName == chromName) {
            chr.blockIds.push_back(blockIdEncoded);
            found = true;
            break;
        }
    }
    if (!found) {
        std::vector< int64_t > ids = { blockIdEncoded };
        tree.chrList.emplace_back((int64_t)chrIdx, ids);
        tree.chrList.back().chrName = chromName;
    }
}

// ----- Chromosome filter parser -----
std::set< std::string > parseChromFilter(const std::string& filterStr) {
    std::set< std::string > result;
    if (filterStr.empty()) return result;
    std::istringstream iss(filterStr);
    std::string token;
    while (std::getline(iss, token, ',')) {
        std::string trimmed;
        for (char c : token) {
            if (!std::isspace(static_cast<unsigned char>(c)))
                trimmed += c;
        }
        if (trimmed.empty()) continue;
        size_t dashPos = trimmed.find('-', 3);
        if (dashPos != std::string::npos && trimmed.rfind("chr", 0) == 0) {
            std::string prefix = "chr";
            try {
                int lo = std::stoi(trimmed.substr(3, dashPos - 3));
                int hi = std::stoi(trimmed.substr(dashPos + 1));
                if (lo > hi) std::swap(lo, hi);
                for (int c = lo; c <= hi; c++) {
                    result.insert(prefix + std::to_string(c));
                }
            } catch (...) {
                result.insert(trimmed);
            }
        } else {
            result.insert(trimmed);
        }
    }
    return result;
}

// ----- Tar handling for chromosome archives -----
// Resolve the per-chromosome MSA archive. Dataset builds use different name
// prefixes (e.g. hprcv2.0_, hprcv1.1_), so accept any file ending in
// "_<chrom>_full_msa.tar.gz" (token-safe: chr1 does not match chr10), a
// bare "<chrom>_full_msa.tar.gz", or the legacy hprcv2.0_ name.
fs::path resolveChromTar(const std::string& tarDir, const std::string& chrom) {
    std::error_code ec;
    fs::path legacy = fs::path(tarDir) / ("hprcv2.0_" + chrom + "_full_msa.tar.gz");
    if (fs::exists(legacy, ec)) return legacy;
    fs::path noPrefix = fs::path(tarDir) / (chrom + "_full_msa.tar.gz");
    if (fs::exists(noPrefix, ec)) return noPrefix;

    const std::string suffix = "_" + chrom + "_full_msa.tar.gz";
    if (fs::exists(fs::path(tarDir), ec)) {
        for (const auto& e : fs::directory_iterator(fs::path(tarDir), ec)) {
            if (!e.is_regular_file(ec)) continue;
            const std::string fn = e.path().filename().string();
            if (fn.size() >= suffix.size() &&
                fn.compare(fn.size() - suffix.size(), suffix.size(), suffix) == 0) {
                return e.path();
            }
        }
    }
    return fs::path();
}

bool extractChromTar(const std::string& tarDir, const std::string& chrom,
                     const std::string& msaRoot) {
    fs::path tarPath = resolveChromTar(tarDir, chrom);
    std::error_code ec;
    if (tarPath.empty() || !fs::exists(tarPath, ec)) {
        std::cerr << "Warning: no MSA tar found for " << chrom << " in " << tarDir
                  << " (expected a file ending in \"_" << chrom
                  << "_full_msa.tar.gz\")" << std::endl;
        return false;
    }

    // Auto-detect strip count by reading first entry
    std::string detectCmd = "tar -tzf " + tarPath.string() + " 2>/dev/null | head -1";
    std::array< char, 512 > buf;
    std::string firstEntry;
    FILE* pipe = popen(detectCmd.c_str(), "r");
    if (pipe) {
        while (fgets(buf.data(), buf.size(), pipe) != nullptr) {
            firstEntry += buf.data();
        }
        pclose(pipe);
    }
    while (!firstEntry.empty() && (firstEntry.back() == '\n' || firstEntry.back() == '\r'))
        firstEntry.pop_back();

    int stripCount = 0;
    if (!firstEntry.empty()) {
        // Find prefix before chrXX/
        std::string chromSlash = chrom + "/";
        size_t chromPos = firstEntry.find(chromSlash);
        if (chromPos != std::string::npos && chromPos > 0) {
            std::string prefix = firstEntry.substr(0, chromPos);
            for (char c : prefix) {
                if (c == '/') stripCount++;
            }
        }
    }

    std::string cmd = "tar -xzf " + tarPath.string();
    if (stripCount > 0) {
        cmd += " --strip-components=" + std::to_string(stripCount);
    }
    cmd += " -C " + msaRoot;

    std::cout << "Extracting " << tarPath.filename().string()
              << " into " << msaRoot;
    if (stripCount > 0) std::cout << " (stripping " << stripCount << " components)";
    std::cout << " ..." << std::endl;

    int ret = std::system(cmd.c_str());
    if (ret != 0) {
        std::cerr << "Error: tar extraction failed (exit " << ret << "): "
                  << cmd << std::endl;
        return false;
    }
    return true;
}

void cleanupExtractedChrom(const std::string& msaRoot, const std::string& chrom) {
    std::error_code ec;
    fs::path msaDir = fs::path(msaRoot) / chrom / "full_msa";
    if (!fs::exists(msaDir, ec)) return;
    std::cout << "Cleaning up extracted .fa files for " << chrom << std::endl;
    size_t removed = 0;
    for (auto& entry : fs::directory_iterator(msaDir, ec)) {
        if (!entry.is_regular_file(ec)) continue;
        std::string ext = entry.path().extension().string();
        if (ext == ".fa" || ext == ".fasta") {
            fs::remove(entry.path(), ec);
            removed++;
        }
    }
    if (removed > 0) {
        std::cout << "  removed " << removed << " extracted MSA files" << std::endl;
    }
}

bool chromMsaDirExists(const std::string& msaRoot, const std::string& chrom) {
    std::error_code ec;
    fs::path nested = fs::path(msaRoot) / chrom / "full_msa";
    fs::path flat   = fs::path(msaRoot) / (chrom + "_full_msa");
    return fs::exists(nested, ec) || fs::exists(flat, ec);
}

} // namespace (end of helpers for buildNetworkFromSprDirs)

panmanUtils::TreeGroup* panmanUtils::buildNetworkFromSprDirs(
    const std::string& treesDir, const std::string& msaRoot,
    const std::string& refFile, const std::string& chromFilter,
    const std::string& tarDir) {

    std::error_code ec;

    // ----- 1. Discover windows from flat Newick files -----
    std::vector< TkWindow > windows;
    if (!fs::exists(treesDir, ec)) {
        std::cerr << "Error: trees directory does not exist: " << treesDir << std::endl;
        return nullptr;
    }
    for (const auto& entry : fs::directory_iterator(treesDir, ec)) {
        if (!entry.is_regular_file(ec)) continue;
        std::string fname = entry.path().filename().string();
        if (!tkEndsWith(fname, ".nwk")) continue;
        std::string baseName = fname.substr(0, fname.size() - 4);
        std::string chrom;
        long long a = 0, b = 0;
        if (!parseWindowName(baseName, chrom, a, b)) continue;
        TkWindow w;
        w.name = baseName;
        w.chrom = chrom;
        w.start = a;
        w.end = b;
        w.newickPath = entry.path().string();
        windows.push_back(w);
    }
    if (windows.empty()) {
        std::cerr << "Error: no .nwk files found in " << treesDir << std::endl;
        return nullptr;
    }

    // ----- 2. Apply chromosome filter -----
    std::set< std::string > chromFilterSet = parseChromFilter(chromFilter);
    if (!chromFilterSet.empty()) {
        std::vector< TkWindow > filtered;
        for (auto& w : windows) {
            if (chromFilterSet.count(w.chrom)) filtered.push_back(w);
        }
        if (filtered.empty()) {
            std::cerr << "Error: no windows match chromosome filter '"
                      << chromFilter << "'" << std::endl;
            return nullptr;
        }
        std::cout << "Chromosome filter: " << windows.size() << " windows -> "
                  << filtered.size() << " windows" << std::endl;
        windows = std::move(filtered);
    }

    // Sort by chromosome then start position
    std::sort(windows.begin(), windows.end(),
              [](const TkWindow& a, const TkWindow& b) {
                  if (a.chrom != b.chrom) return a.chrom < b.chrom;
                  return a.start < b.start;
              });

    // Assign chromosome indices
    std::map< std::string, int > chromToIdx;
    for (auto& w : windows) {
        if (chromToIdx.find(w.chrom) == chromToIdx.end()) {
            int idx = static_cast<int>(chromToIdx.size());
            chromToIdx[w.chrom] = idx;
        }
        w.chrIdx = chromToIdx[w.chrom];
    }

    std::cout << "Discovered " << windows.size() << " windows across "
              << chromToIdx.size() << " chromosomes." << std::endl;

    // ----- 3. Load reference chromosomes -----
    std::unordered_map< std::string, std::string > refChromCache;
    auto getRefSeq = [&](const std::string& chrom) -> const std::string& {
        auto it = refChromCache.find(chrom);
        if (it != refChromCache.end()) return it->second;
        std::string seq = refFile.empty() ? std::string()
                                          : loadReferenceChromosome(refFile, chrom);
        refChromCache[chrom] = std::move(seq);
        return refChromCache[chrom];
    };

    fs::path tmpDir = fs::path(".") / "tmp_spr_build";
    fs::create_directories(tmpDir, ec);

    // ----- 4. Infer SPR moves between consecutive windows -----
    std::vector< SprStep > sprSteps = collectSprSteps(windows);
    std::cout << "Inferred " << sprSteps.size() << " SPR moves from "
              << windows.size() << " windows." << std::endl;

    // ----- 5. Build single tree from window 0 -----
    std::string extractedChrom;
    auto ensureChromExtracted = [&](const std::string& chrom) {
        if (tarDir.empty()) return;
        if (extractedChrom == chrom) return;
        if (!extractedChrom.empty()) {
            cleanupExtractedChrom(msaRoot, extractedChrom);
            extractedChrom.clear();
        }
        if (extractChromTar(tarDir, chrom, msaRoot)) {
            extractedChrom = chrom;
        }
    };
    auto cleanupCurrentChrom = [&]() {
        if (!extractedChrom.empty()) {
            cleanupExtractedChrom(msaRoot, extractedChrom);
            extractedChrom.clear();
        }
    };

    const TkWindow& anchorWin = windows[0];
    ensureChromExtracted(anchorWin.chrom);

    std::cout << "Building base tree from " << anchorWin.name
              << " (window 1/" << windows.size() << ")" << std::endl;
    panmanUtils::Tree* tree0Ptr =
        buildWindowPanmat(anchorWin, windows, 0, msaRoot, tmpDir, getRefSeq);
    if (!tree0Ptr) {
        std::cerr << "Error: failed to build base tree from " << anchorWin.name
                  << std::endl;
        cleanupCurrentChrom();
        fs::remove_all(tmpDir, ec);
        return nullptr;
    }
    annotateBlockChromosome(*tree0Ptr, anchorWin.chrIdx, anchorWin.chrom);
    if (tree0Ptr->root) {
        tree0Ptr->root->isComMutHead = true;
        tree0Ptr->root->treeIndex = 0;
    }
    if (sprDebugEnabled()) {
        sprDebugTreeStats("base tree (window 0)", *tree0Ptr);
    }

    // Build a lookup: windowIdx → SprStep
    std::unordered_map< size_t, size_t > windowToSprStep;
    for (size_t s = 0; s < sprSteps.size(); s++) {
        windowToSprStep[sprSteps[s].windowIdx] = s;
    }

    // ----- 6. Process windows 1..N -----
    std::vector< panmanUtils::SerializedNetworkEdge > networkEdges;
    std::vector< AppliedSpr > appliedSprs;
    size_t sprBuilt = 0;
    size_t sprSkippedMap = 0;
    size_t sprSkippedCycle = 0;
    size_t blocksAdded = 1;

    // Union graph = base-tree edges (parent -> child) plus every emitted network
    // edge (acceptor -> moved). A tskit tree sequence needs one global time per
    // node with time[parent] > time[child] on every edge, so this union must be
    // acyclic. We seed it with the base tree and reject any SPR network edge that
    // would close a directed cycle (keeping the earlier edge as the parent).
    std::unordered_map< std::string, std::vector< std::string > > unionChildren;
    for (const auto& kv : tree0Ptr->allNodes) {
        const panmanUtils::Node* n = kv.second;
        if (n->parent) unionChildren[n->parent->identifier].push_back(n->identifier);
    }
    // Would adding edge acceptor -> moved create a cycle? It does iff `moved`
    // can already reach `acceptor` through the current union graph.
    auto wouldCloseCycle = [&](const std::string& moved,
                               const std::string& acceptor) -> bool {
        if (moved == acceptor) return true;
        std::unordered_set< std::string > visited;
        std::vector< std::string > stack{ moved };
        visited.insert(moved);
        while (!stack.empty()) {
            std::string cur = std::move(stack.back());
            stack.pop_back();
            auto it = unionChildren.find(cur);
            if (it == unionChildren.end()) continue;
            for (const std::string& nxt : it->second) {
                if (nxt == acceptor) return true;
                if (visited.insert(nxt).second) stack.push_back(nxt);
            }
        }
        return false;
    };

    auto removeUnionChild = [&](const std::string& parent,
                                const std::string& child) {
        auto it = unionChildren.find(parent);
        if (it == unionChildren.end()) return;
        auto& ch = it->second;
        ch.erase(std::remove(ch.begin(), ch.end(), child), ch.end());
    };

    std::mt19937 cycleRng(std::random_device{}());

    for (size_t wIdx = 1; wIdx < windows.size(); wIdx++) {
        const TkWindow& wCurr = windows[wIdx];

        ensureChromExtracted(wCurr.chrom);

        // 6a. Apply SPR if one was inferred
        auto sprIt = windowToSprStep.find(wIdx);
        if (sprIt != windowToSprStep.end()) {
            const SprStep& step = sprSteps[sprIt->second];
            const SprMoveInfo& spr = step.spr;
            const std::string edgeLabel =
                windows[wIdx - 1].name + " -> " + wCurr.name;

            if (sprDebugEnabled()) {
                std::ostringstream os;
                os << "--- SPR at window " << wIdx << " " << edgeLabel;
                sprDebugLog(os.str());
            }

            panmanUtils::Node* movedNode =
                findMRCA(*tree0Ptr, spr.movedClade);
            panmanUtils::Node* acceptorNode =
                findMRCA(*tree0Ptr, spr.acceptorClade);
            panmanUtils::Node* donorNode =
                findMRCA(*tree0Ptr, spr.donorClade);

            if (!acceptorNode && donorNode && donorNode->parent) {
                acceptorNode = donorNode->parent;
            }

            bool wouldCycleInTree = movedNode && acceptorNode &&
                                    isAncestorOf(movedNode, acceptorNode);
            if (!movedNode || !acceptorNode ||
                movedNode == acceptorNode || movedNode == tree0Ptr->root ||
                wouldCycleInTree) {
                if (sprDebugEnabled()) {
                    sprDebugLog(wouldCycleInTree
                        ? "skip SPR: acceptor is descendant of moved in current tree"
                        : "skip SPR: mapping failed or trivial");
                    sprDebugClade("moved", spr.movedClade);
                    sprDebugClade("acceptor", spr.acceptorClade);
                }
                sprSkippedMap++;
            } else {
                const std::string& accId = acceptorNode->identifier;
                const std::string& movId = movedNode->identifier;
                bool emitEdge = true;

                if (wouldCloseCycle(movId, accId)) {
                    // Direct 2-cycle: union already has movId -> accId and this SPR
                    // wants accId -> movId. Randomly keep one parent direction.
                    bool hasReverse = false;
                    auto revIt = unionChildren.find(movId);
                    if (revIt != unionChildren.end()) {
                        for (const std::string& c : revIt->second) {
                            if (c == accId) { hasReverse = true; break; }
                        }
                    }
                    if (hasReverse && (cycleRng() & 1)) {
                        networkEdges.erase(
                            std::remove_if(networkEdges.begin(), networkEdges.end(),
                                [&](const panmanUtils::SerializedNetworkEdge& e) {
                                    return e.parentNodeId == movId &&
                                           e.childNodeId == accId;
                                }),
                            networkEdges.end());
                        removeUnionChild(movId, accId);
                        if (sprDebugEnabled()) {
                            sprDebugLog("cycle resolve: kept new parent " + accId +
                                        " -> " + movId + " (dropped reverse)");
                        }
                    } else {
                        emitEdge = false;
                        sprSkippedCycle++;
                        if (sprDebugEnabled()) {
                            sprDebugLog(hasReverse
                                ? "cycle resolve: kept existing parent " + movId +
                                  " -> " + accId + " (skipped new SPR edge)"
                                : "skip SPR: would close directed cycle in union graph");
                        }
                    }
                }

                if (emitEdge) {
                    networkEdges.emplace_back(
                        accId, movId,
                        std::vector< int64_t >{ static_cast<int64_t>(wIdx) });
                    unionChildren[accId].push_back(movId);

                    appliedSprs.push_back(
                        applySprToTree(*tree0Ptr, movedNode, acceptorNode));

                    sprBuilt++;
                    if (sprDebugEnabled()) {
                        std::ostringstream os;
                        os << "applied SPR: reparent "
                           << movId << " -> " << accId;
                        sprDebugLog(os.str());
                    }
                }
            }
        }

        // 6b. Read MSA for this window
        std::string msaPath = resolveMsaPath(msaRoot, wCurr);
        if (msaPath.empty()) {
            std::cerr << "Warning: no MSA found for " << wCurr.name
                      << "; skipping block." << std::endl;
            continue;
        }

        std::string useMsaPath = msaPath;
        std::string tmpMsa;
        std::string suffix;
        if (wIdx + 1 < windows.size() &&
            windows[wIdx + 1].chrom == wCurr.chrom) {
            const std::string& refSeq = getRefSeq(wCurr.chrom);
            if (!refSeq.empty()) {
                long long gapBeg = wCurr.end;
                long long gapLen = windows[wIdx + 1].start - wCurr.end - 1;
                if (gapLen > 0 && gapBeg >= 0 &&
                    gapBeg + gapLen <=
                        static_cast<long long>(refSeq.size())) {
                    suffix = refSeq.substr(static_cast<size_t>(gapBeg),
                                           static_cast<size_t>(gapLen));
                }
            }
        }
        if (!suffix.empty()) {
            tmpMsa = (tmpDir / (wCurr.name + ".aug.fa")).string();
            if (writeAugmentedMSA(msaPath, suffix, tmpMsa)) {
                useMsaPath = tmpMsa;
            } else {
                tmpMsa.clear();
            }
        }

        std::map< std::string, std::string > leafSeqs;
        size_t alignLen = readMsaFile(useMsaPath, leafSeqs);
        if (!tmpMsa.empty()) fs::remove(tmpMsa, ec);

        if (leafSeqs.empty() || alignLen == 0) {
            std::cerr << "Warning: empty MSA for " << wCurr.name
                      << "; skipping block." << std::endl;
            continue;
        }

        // 6c. Run Fitch on the current (SPR-modified) topology
        int32_t blockId = static_cast<int32_t>(wIdx);
        runFitchForBlock(*tree0Ptr, leafSeqs, alignLen, blockId,
                         std::string("GRCh38"));
        annotateBlockChromosome(*tree0Ptr, wCurr.chrIdx, wCurr.chrom,
                                blockId);
        blocksAdded++;

        std::cout << "Window " << (wIdx + 1) << "/" << windows.size()
                  << " " << wCurr.name
                  << ": block " << blockId
                  << " (" << alignLen << " cols, "
                  << leafSeqs.size() << " seqs)" << std::endl;
    }

    // Cleanup any remaining extracted chromosome
    cleanupCurrentChrom();

    // ----- 7. Restore base tree topology -----
    for (auto it = appliedSprs.rbegin(); it != appliedSprs.rend(); ++it) {
        undoSpr(*it);
    }
    appliedSprs.clear();

    if (sprDebugEnabled()) {
        std::ostringstream os;
        os << "SPR summary: applied=" << sprBuilt
           << " skipped=" << sprSkippedMap
           << " cycle_skipped=" << sprSkippedCycle
           << " edges=" << networkEdges.size()
           << " blocks=" << blocksAdded;
        sprDebugLog(os.str());
    }

    fs::remove_all(tmpDir, ec);

    // ----- 8. Create TreeGroup with 1 tree + network edges -----
    std::vector< panmanUtils::Tree* > treePtrs = { tree0Ptr };
    panmanUtils::TreeGroup* TG = new panmanUtils::TreeGroup(treePtrs);
    TG->networkEdges = std::move(networkEdges);

    std::cout << "Built single-tree network: " << windows.size() << " windows, "
              << blocksAdded << " blocks, "
              << tree0Ptr->allNodes.size() << " nodes, "
              << chromToIdx.size() << " chromosomes, "
              << TG->networkEdges.size() << " network edges";
    if (sprSkippedCycle > 0) {
        std::cout << " (" << sprSkippedCycle
                  << " SPR edge(s) skipped to keep union graph acyclic)";
    }
    std::cout << "." << std::endl;

    delete tree0Ptr;
    return TG;
}
