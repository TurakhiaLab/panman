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

#include <functional>
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
    std::string newickPath;  // path to <name>_resolved.nwk
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
