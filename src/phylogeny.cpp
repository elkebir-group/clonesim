//
// Created by Mohammed El-Kebir on 8/20/23.
//

#include "phylogeny.h"
#include <fstream>
#include <algorithm>
#include <lemon/connectivity.h>
#include "beta_distribution.hpp"
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/discrete_distribution.hpp>

Phylogeny::Phylogeny()
        : _T(), _root(_T.addNode()), _mrca(lemon::INVALID), _trunk(_T), _cnaTrees(), _charState(_T), _trunkLength(0),
          _mutToCluster(), _clusterToMut(), _xbar(_T), _ybar(_T), _clusterToNode(), _nodeToCluster(_T), _mutToSegment(),
          _segmentToMut(), _D(_T), _clusterD(), _proportions(_T), _nodeToIndex(_T) {
    _trunk[_root] = false;
}

Phylogeny::Phylogeny(const Phylogeny &other)
        : _T(), _root(lemon::INVALID), _mrca(lemon::INVALID), _trunk(_T), _cnaTrees(other._cnaTrees), _charState(_T),
          _trunkLength(other._trunkLength), _mutToCluster(other._mutToCluster), _clusterToMut(other._clusterToMut),
          _xbar(_T), _ybar(_T), _clusterToNode(other._clusterToNode), _nodeToCluster(_T),
          _mutToSegment(other._mutToSegment), _segmentToMut(other._segmentToMut), _D(_T), _clusterD(), _proportions(_T),
          _nodeToIndex(_T) {
    lemon::digraphCopy(other._T, _T)
            .node(other._root, _root)
            .node(other._mrca, _mrca)
            .nodeMap(other._trunk, _trunk)
            .nodeMap(other._charState, _charState)
            .nodeMap(other._xbar, _xbar)
            .nodeMap(other._ybar, _ybar)
            .nodeMap(other._nodeToCluster, _nodeToCluster)
            .nodeMap(other._proportions, _proportions)
            .nodeMap(other._nodeToIndex, _nodeToIndex)
            .nodeMap(other._D, _D)
                    //.nodeMatrix(other._clusterD, _clusterD)

                    //.arcMap(other._T.SubArcIt)
            .run();

    _clusterToNode = NodeVector(_clusterToMut.size(), lemon::INVALID);

    for (NodeIt v(_T); v != lemon::INVALID; ++v) {
        if (_nodeToCluster[v] != -1) {
            _clusterToNode[_nodeToCluster[v]] = v;
        }
    }


    initD(_root);
    _clusterD = NodeMatrix(_clusterToNode.size(), NodeVector(0));
    //initClusterD(_mrca, 0);
    initClusterD();
}

void Phylogeny::addSegment(const CnaTree &cnaTree, bool truncal) {
    _cnaTrees.push_back(cnaTree);
    _cnaTrees.back() = cnaTree;
    _segmentToMut.push_back(IntSet());

    Digraph G;
    Node rootG = lemon::INVALID;
    NodePairNodeMap mapG(G);
    NodeNodeSetMap toG(_T);
    NodeNodeSetMap cnaToG(cnaTree.T());

    createProductGraph(cnaTree, G, mapG, toG, cnaToG, rootG);

    BoolNodeMap filterNodesG(G, true);
    BoolArcMap filterArcsG(G, true);
    SubDigraph subG(G, filterNodesG, filterArcsG);
//  printProductGraph(G, subG, mapG, cnaTree, std::cout);

    sample(G, mapG, toG, cnaToG, rootG, subG);

//  printProductGraph(G, subG, mapG, cnaTree, std::cerr);

    update(cnaTree, truncal, G, mapG, toG, cnaToG, rootG, subG);

//  writeDOT(std::cerr);
//  writeDOT(std::cerr);

//  assert(lemon::countStronglyConnectedComponents(_T) == 1);

    initD(_root);

    assert(lemon::dag(_T));
}

void Phylogeny::update(const CnaTree &cnaTree, bool truncal,
                       const Digraph &G, const NodePairNodeMap &mapG, const NodeNodeSetMap &toG,
                       const NodeNodeSetMap &cnaToG, const Node rootG, const SubDigraph &subG) {
    NodeSet toErase;
    for (NodeIt vv(_T); vv != lemon::INVALID; ++vv) {
        toErase.insert(vv);
    }

    NodeNodeMap productNodesToNewNodes(G, lemon::INVALID);
    for (SubNodeIt uv(subG); uv != lemon::INVALID; ++uv) {
        if (SubOutArcIt(subG, uv) != lemon::INVALID || SubInArcIt(subG, uv) != lemon::INVALID) {
            Node new_uv = _T.addNode();
            _nodeToCluster[new_uv] = -1;
            _charState[new_uv] = _charState[mapG[uv].first];
            _charState[new_uv].push_back(cnaTree.state(mapG[uv].second));
            _trunk[new_uv] = false;
//      _trunk[new_uv] = _trunk[mapG[uv].first];
            productNodesToNewNodes[uv] = new_uv;

            if (mapG[uv].first == _root && cnaTree.state(mapG[uv].second) == cnaTree.rootState()) {
                _root = new_uv;
            }
//      if ((mapG[uv].first == _mrca || (_mrca == lemon::INVALID && mapG[uv].first == _root)) && cnaTree.truncal() && truncal)
////
////        && mapG[uv].second == cnaTree.T().target(OutArcIt(cnaTree.T(), mapG[uv].second)))
//      {
//        _mrca = new_uv;
//      }

        }
    }

//
//  // kind of a hack
//  _trunk[_root] = true;

    for (SubArcIt aa(subG); aa != lemon::INVALID; ++aa) {
        Node source = subG.source(aa);
        Node target = subG.target(aa);
        Node new_source = productNodesToNewNodes[source];
        Node new_target = productNodesToNewNodes[target];
        _T.addArc(new_source, new_target);
//    if (new_target == _mrca)
//    {
//      _trunk[new_source] = true;
//    }
    }

    // determine new trunk and mrca
    _mrca = _root;
    int trunkLength = 0;
    while (lemon::countOutArcs(_T, _mrca) == 1) {
        _trunk[_mrca] = true;
        _mrca = _T.target(OutArcIt(_T, _mrca));
        ++trunkLength;
        if (trunkLength == _trunkLength + 1)
            break;
    }
    _trunkLength = trunkLength;

    for (Node node: toErase) {
        _T.erase(node);
    }
}

void Phylogeny::createIndex() {
    // Digraph::NodeMap<int> nodeToIndex(_T);
    int idx = 0;
    for (NodeIt v(_T); v != lemon::INVALID; ++v) {
        _nodeToIndex[v] = idx++;
    }
    // Iterate over the nodes in the map and print their associated values
    // for ( NodeIt node(_T); node != lemon::INVALID; ++node) {
    //    int value = _nodeToIndex[node]; // Access the integer value associated with the node
    //   std::cout << "Node: " << value << std::endl;
    //}

}

void Phylogeny::sample(const Digraph &G, const NodePairNodeMap &mapG, const NodeNodeSetMap &toG,
                       const NodeNodeSetMap &cnaToG, const Node rootG, SubDigraph &subG) const {
    // build frontier
    ArcSet F;
    for (OutArcIt a(G, rootG); a != lemon::INVALID; ++a) {
        F.insert(a);
    }

    ArcSet sampledEdges;
    while (!F.empty()) {
        // pick edge at random
        Arc a = lemon::INVALID;

        int sample = boost::random::uniform_int_distribution<>(1, F.size())(g_rng);

        int i = 1;
        for (Arc aa: F) {
            if (i++ == sample) {
                a = aa;
                break;
            }
        }
        assert(a != lemon::INVALID);
        F.erase(a);
        sampledEdges.insert(a);

        Node uv_source = G.source(a);
        Node uv_target = G.target(a);

        assert((mapG[uv_source].first == mapG[uv_target].first) || (mapG[uv_source].second == mapG[uv_target].second));

        if (mapG[uv_source].first == mapG[uv_target].first) {
            // remove all incoming arcs to mapG[uv_target].second
            for (Node uv: cnaToG[mapG[uv_target].second]) {
                for (InArcIt in_uv(G, uv); in_uv != lemon::INVALID; ++in_uv) {
                    if (sampledEdges.count(in_uv) == 0 && (mapG[G.source(in_uv)].second == mapG[uv_source].second) &&
                        (mapG[G.target(in_uv)].second == mapG[uv_target].second)) {
                        subG.disable(in_uv);
                        F.erase(in_uv);
                    }
                }
            }
        }

        if (mapG[uv_source].second == mapG[uv_target].second) {
            // remove all incoming arcs to mapG[uv_target].first
            for (Node uv: toG[mapG[uv_target].first]) {
                for (InArcIt in_uv(G, uv); in_uv != lemon::INVALID; ++in_uv) {
                    if (sampledEdges.count(in_uv) == 0 && (mapG[G.source(in_uv)].first == mapG[uv_source].first) &&
                        (mapG[G.target(in_uv)].first == mapG[uv_target].first)) {
                        subG.disable(in_uv);
                        F.erase(in_uv);
                    }
                }
            }
        }

        for (SubOutArcIt aa(subG, uv_target); aa != lemon::INVALID; ++aa) {
            F.insert(aa);
        }
    }
}

void Phylogeny::createProductGraph(const CnaTree &cnaTree, Digraph &G, NodePairNodeMap &mapG,
                                   NodeNodeSetMap &toG, NodeNodeSetMap &cnaToG, Node &rootG) const {
    // add nodes to product graph
    for (NodeIt u(_T); u != lemon::INVALID; ++u) {
        for (NodeIt v(cnaTree.T()); v != lemon::INVALID; ++v) {
            Node uv = G.addNode();
            mapG[uv] = std::make_pair(u, v);
            toG[u].insert(uv);
            cnaToG[v].insert(uv);

            if (u == _root && cnaTree.state(v) == cnaTree.rootState()) {
                rootG = uv;
            }
        }
    }

    // add edges to product graph
    for (ArcIt a(_T); a != lemon::INVALID; ++a) {
        for (Node uv: toG[_T.source(a)]) {
            for (Node uv2: toG[_T.target(a)]) {
                if (mapG[uv].second == mapG[uv2].second) {
                    G.addArc(uv, uv2);
                }
            }
        }
    }
    for (ArcIt a(cnaTree.T()); a != lemon::INVALID; ++a) {
        for (Node uv: cnaToG[cnaTree.T().source(a)]) {
            for (Node uv2: cnaToG[cnaTree.T().target(a)]) {
                if (mapG[uv].first == mapG[uv2].first && !_trunk[mapG[uv].first]) {
                    G.addArc(uv, uv2);
                }
            }
        }
    }
}

void Phylogeny::printProductGraph(const Digraph &G,
                                  const SubDigraph &subG, const NodePairNodeMap &mapG, const CnaTree &cnaTree,
                                  std::ostream &out) const {
    out << "digraph G {" << std::endl;

    for (SubNodeIt v(subG); v != lemon::INVALID; ++v) {
        out << "\t" << subG.id(v) << " [label=\"";
        int idx = 0;
        for (int state: _charState[mapG[v].first]) {
            out << _cnaTrees[idx].label(state) << ",";
            idx++;
        }
        out << cnaTree.label(cnaTree.state(mapG[v].second));
        if (mapG[v].first == _root) {
            out << "\\nroot";
        }
        if (mapG[v].first == _mrca) {
            out << "\\nmrca";
        }
        if (_trunk[mapG[v].first]) {
            out << "\\ntrunk";
        }
        out << "\"";
        if (!(SubOutArcIt(subG, v) != lemon::INVALID || SubInArcIt(subG, v) != lemon::INVALID)) {
            out << ",style=dashed";
        }
        if (_trunk[mapG[v].first]) {
            out << ",color=red";
        }
        out << "]" << std::endl;
    }

    for (ArcIt a(G); a != lemon::INVALID; ++a) {
        out << "\t" << G.id(G.source(a)) << " -> " << G.id(G.target(a));
        if (!subG.status(a)) {
            out << " [style=dashed]";
        }
        out << std::endl;
    }

    out << "}" << std::endl;
}

void Phylogeny::writeDOT(std::ostream &out) const {
    out << "digraph T {" << std::endl;
    for (NodeIt v(_T); v != lemon::INVALID; ++v) {
        out << "\t" << _T.id(v) << " [label=\"";
        for (int i = 0; i < _charState[v].size(); ++i) {
//      out << _cnaTrees[i].label(_charState[v][i]) << " ; ";
            out << _cnaTrees[i].label(_charState[v][i]) << " ; ";
        }
        for (int segmentIdx = 0; segmentIdx < _cnaTrees.size(); ++segmentIdx) {
            out << "\\n";
            for (int mutIdx: _segmentToMut[segmentIdx]) {
                out << _xbar[v][mutIdx] << "|" << _ybar[v][mutIdx] << " ; ";
            }
        }
        if (v == _root) {
            out << "\\nroot";
        }
        if (v == _mrca) {
            out << "\\nmrca";
        }
        if (_trunk[v]) {
            out << "\\ntrunk";
        }
        if (_nodeToCluster[v] != -1) {
            out << "\\ncluster " << _nodeToCluster[v];
        }
        out << "\\n[";
        for (double prop: _proportions[v]) {
            out << " " << prop;
        }
        out << " ]";
        out << "\"]" << std::endl;
    }

    for (ArcIt a(_T); a != lemon::INVALID; ++a) {
        out << _T.id(_T.source(a)) << " -> " << _T.id(_T.target(a)) << std::endl;
    }

    out << "}" << std::endl;
}

void Phylogeny::writeFiles(std::ostream &out, std::string &outputDir, int nrsamples) const {
    std::string outputTreeFile = outputDir + "/tree.tsv";
    std::string outputNodeFile = outputDir + "/node.tsv";
    std::string outputProportionFile = outputDir + "/proportions.tsv";
    writeTree(outputTreeFile);
    writeNodeFile(out, outputNodeFile);
    writeProportionFile(out, outputProportionFile, nrsamples);
}

void Phylogeny::writeTree(std::string &treefn) const {
    Digraph::NodeMap<int> nodeToIndex_updated(_T);
    int idx = 0;
    for (NodeIt v(_T); v != lemon::INVALID; ++v) {
        nodeToIndex_updated[v] = idx++;
    }
    std::ofstream myfile;
    if (treefn.compare("") != 0) {
        myfile.open(treefn);
    }
    if (!myfile) {
        throw std::runtime_error("Error in creating output file");
    }
    for (ArcIt a(_T); a != lemon::INVALID; ++a) {
        Node par = _T.source(a);
        int parent = nodeToIndex_updated[par];
        Node ch = _T.target(a);
        int child = nodeToIndex_updated[ch];
        myfile << parent << "\t" << child << std::endl;
    }
    myfile.close();

}

void Phylogeny::writeNodeFile(std::ostream &out, std::string &outputNodeFilename) const {
    Digraph::NodeMap<int> nodeToIndex_updated(_T);
    int idx = 0;
    for (NodeIt v(_T); v != lemon::INVALID; ++v) {
        nodeToIndex_updated[v] = idx++;
    }
    if (outputNodeFilename.compare("") != 0) {
        std::ofstream myfile;
        myfile.open(outputNodeFilename);
        if (!myfile) {
            throw std::runtime_error("Error in creating output file");
        } else {

            int isThereMutation = 0;
            myfile << "node\tsegment\tx\ty\tm\txbar\tybar\tclusterIDofMutation\n";
            for (NodeIt v(_T); v != lemon::INVALID; ++v) // printing out node
            {
                for (int segmentIdx = 0; segmentIdx < _cnaTrees.size(); ++segmentIdx) {
                    for (int mutIdx: _segmentToMut[segmentIdx]) {
                        isThereMutation = 1;
                        myfile << nodeToIndex_updated[v] << "\t"; //printing out node
                        myfile << segmentIdx << "\t"; //segment label
                        myfile << _cnaTrees[segmentIdx].label(_charState[v][segmentIdx])[1] << "\t"; //x of segment
                        myfile << _cnaTrees[segmentIdx].label(_charState[v][segmentIdx])[3] << "\t"; //y of segment
                        myfile << mutIdx << "\t"; //mutationID
                        myfile << _xbar[v][mutIdx] << "\t" << _ybar[v][mutIdx] << "\t"; //xbar and ybar of mutation
                        myfile << _mutToCluster[mutIdx] << "\n"; //clusterID of mutation
                    }
                    if (isThereMutation == 0) { //edge case for if no mutations occur in that segment
                        myfile << nodeToIndex_updated[v] << "\t"; //printing out node
                        myfile << segmentIdx << "\t"; //segment label
                        myfile << _cnaTrees[segmentIdx].label(_charState[v][segmentIdx])[1] << "\t"; //x of segment
                        myfile << _cnaTrees[segmentIdx].label(_charState[v][segmentIdx])[3] << "\t"; //y of segment
                        myfile << -1 << "\t"; //mutationID
                        myfile << -1 << "\t" << -1 << "\t"; //xbar and ybar of mutation
                        myfile << -1 << "\n"; //clusterID of mutation
                    }
                    isThereMutation = 0;
                }
            }

        }
        myfile.close();
    } else {
        throw std::runtime_error("Please specify proportion file name");
    }


}

std::vector<float> dirichlet(double dirchparam, int n, bool uniform, float threshold)
{

    if (!uniform)
    {

        boost::random::gamma_distribution<> gamma_dist(dirchparam, 1); //assuming scaleparam of 1
        std::vector<float> gamma_sample(n);

        //sample from the dirichelet distribution until all values are above the threshold
        bool accept = true;
        do{
            
            accept = true;
            float gamma_sum = 0;
            for (int i = 0; i < n; i++)
            {

                gamma_sample[i] = gamma_dist(g_rng);

                gamma_sum += gamma_sample[i];
            }

            for (int i = 0; i < n; i++)
            {
                gamma_sample[i] = gamma_sample[i] / gamma_sum;
                if(gamma_sample[i] < threshold){
                    accept = false;
                }
      
            }

            //}
        }while(!accept);


   
        return gamma_sample;

    }
    else
    {
        std::cout <<"uniform sizes"  << std::endl;
        float fraction = 1.0 / n;
        std::vector<float> unf(n, fraction);
        return unf;
    }

}

void Phylogeny::sampleMutations(int n, int l, int num_tries, double dirich_param, bool uniform, double threshold) {
    int validSampling = 0;
    int numberOfTries = 0;
    while (validSampling == 0) {
        numberOfTries += 1;
        if (!(n >=1 & l >= 1 & l <= n)) {
            throw std::runtime_error("Invalid input: number of mutation clusters must be <= number of SNVs, and both must be >= 1");
        }
        if (numberOfTries > num_tries) {
            throw std::runtime_error("Rejection sampling has exceeded the number of possible tries for sampling mutations. Increase the number of tries or try a different seed.");
        }
        validSampling = 1;

        const int k = getNrSegments();
        _clusterToMut = IntSetVector(l);
        _mutToSegment = IntVector();
        _segmentToMut = IntSetVector();
        for (int i = 0; i < k; i++) {
            _segmentToMut.push_back(IntSet());
        }
        _mutToCluster = IntVector();
        _clusterToNode = NodeVector();

        // assign mutations to segments
        //std::uniform_int_distribution<> uniform_segments(0, k - 1);
        boost::random::uniform_01<double> myrand;


        //sample from a weighted discrete distribution to assign each mutation i to one of k segments
        std::vector<float> dir_segments = dirichlet(dirich_param, k, uniform, threshold);
        boost::random::discrete_distribution<int> segment_dist(dir_segments.begin(), dir_segments.end());

        //sample from a weighted discrete distribution to assign each mutation i to one of l clusters
        std::vector<float> dir_clusters = dirichlet(dirich_param, l, uniform, threshold);
        
        //sort the weights in descending order
        std::sort(dir_clusters.rbegin(), dir_clusters.rend());
        boost::random::discrete_distribution<int> cluster_dist(dir_clusters.begin(), dir_clusters.end());
        
        for(int i=0; i < n; ++i)
        {
            int segmentIdx = segment_dist(g_rng);
            _mutToSegment.push_back(segmentIdx);
            _segmentToMut[segmentIdx].insert(i);

            int clusterIdx = cluster_dist(g_rng);
            _mutToCluster.push_back(clusterIdx);
            _clusterToMut[clusterIdx].insert(i);
        }


        NodeVector trunkNodes;
        trunkNodes.push_back(_mrca);
        for (NodeIt v(_T); v != lemon::INVALID; ++v) {
            if (_trunk[v] && v != _root) {
                trunkNodes.push_back(v);
            }
        }

        if (trunkNodes.empty()) {
            trunkNodes.push_back(_root);
        }

        // pick trunk mutation cluster
        boost::random::uniform_int_distribution<> unif1(0, trunkNodes.size() - 1);
        int idx = unif1(g_rng);
        _clusterToNode.push_back(trunkNodes[idx]);
        _nodeToCluster[trunkNodes[idx]] = 0;

        // sample the remaining l-1 mutation cluster locations
        NodeVector remainingNodes(_D[_clusterToNode[0]].begin(), _D[_clusterToNode[0]].end());
        std::shuffle(remainingNodes.begin(), remainingNodes.end(), g_rng);
        if (!(remainingNodes.size() >= l - 1)) {
            throw std::runtime_error("Invalid number of clusters must be smaller than the number of nodes.");
        }

        for (int i = 0; i < l - 1; ++i) {
            _clusterToNode.push_back(remainingNodes[i]);
            _nodeToCluster[remainingNodes[i]] = i + 1;
        }

        boost::random::uniform_01<double> unif01;

        for (NodeIt v(_T); v != lemon::INVALID; ++v) {
            _xbar[v] = IntVector(n, 0);
            _ybar[v] = IntVector(n, 0);
        }

        // sample xbar and ybar for each mutation
        for (int mutIdx = 0; mutIdx < n; ++mutIdx) {
            int clusterIdx = _mutToCluster[mutIdx];
            int segmentIdx = _mutToSegment[mutIdx];
            Node mutationNode = _clusterToNode[clusterIdx];

            int x = _cnaTrees[segmentIdx].x(_charState[mutationNode][segmentIdx]);
            int y = _cnaTrees[segmentIdx].y(_charState[mutationNode][segmentIdx]);

            if (!(x > 0 | y > 0)) {
                validSampling = 0;
            }

            if (x > 0 && y > 0) {
                bool mut_x = unif01(g_rng) > 0.5;
                if (mut_x) {
                    _xbar[mutationNode][mutIdx] = 1;
                    _ybar[mutationNode][mutIdx] = 0;
                    sampleMutation(mutationNode, segmentIdx, mutIdx);
                } else {
                    _xbar[mutationNode][mutIdx] = 0;
                    _ybar[mutationNode][mutIdx] = 1;
                    sampleMutation(mutationNode, segmentIdx, mutIdx);
                }
            } else if (x > 0) {
                _xbar[mutationNode][mutIdx] = 1;
                _ybar[mutationNode][mutIdx] = 0;
                sampleMutation(mutationNode, segmentIdx, mutIdx);
            } else if (y > 0) {
                _xbar[mutationNode][mutIdx] = 0;
                _ybar[mutationNode][mutIdx] = 1;
                sampleMutation(mutationNode, segmentIdx, mutIdx);
            }
        }

    }
    _clusterD = NodeMatrix(_clusterToMut.size(), NodeVector(0));
    initClusterD();


}


void Phylogeny::sampleMutation(const Node u, const int segmentIdx, const int mutIdx) {
    int x_u = _cnaTrees[segmentIdx].x(_charState[u][segmentIdx]);
    int y_u = _cnaTrees[segmentIdx].y(_charState[u][segmentIdx]);
    int xbar_u = _xbar[u][mutIdx];
    int ybar_u = _ybar[u][mutIdx];
    if (! (0 <= xbar_u && xbar_u <= x_u)) {
        throw std::runtime_error("Error in sampleMutation: the variant copy number must be less than the segment copy number for x");
    }
    if (! (0 <= ybar_u && ybar_u <= y_u)) {
        throw std::runtime_error("Error in sampleMutation: the variant copy number must be less than the segment copy number for y");
    }

    for (OutArcIt a(_T, u); a != lemon::INVALID; ++a) {
        Node v = _T.target(a);
        int x_v = _cnaTrees[segmentIdx].x(_charState[v][segmentIdx]);
        int y_v = _cnaTrees[segmentIdx].y(_charState[v][segmentIdx]);

        if (xbar_u == 0 && ybar_u == 0) {
            _xbar[v][mutIdx] = _ybar[v][mutIdx];
        } else if (xbar_u > 0) {
            if (x_u == x_v) {
                _xbar[v][mutIdx] = xbar_u;
                _ybar[v][mutIdx] = ybar_u;
            } else if (xbar_u == x_u) {
                _xbar[v][mutIdx] = x_v;
                _ybar[v][mutIdx] = 0;
            } else if (x_u < x_v) {
                // amplification
                int delta = x_v - x_u;
                IntVector choices;
                for (int l = 0; l <= delta; ++l) {
                    choices.push_back(xbar_u + delta);
                }

                boost::random::uniform_int_distribution<> unif(0, choices.size() - 1);
                //std::uniform_int_distribution<> unif(0, choices.size() - 1);
                _xbar[v][mutIdx] = choices[unif(g_rng)];
                _ybar[v][mutIdx] = 0;
            } else {
                // deletion
                if (!(x_u > x_v)) {
                    throw std::runtime_error("Error in samplemutation: x_u must be > x_v");
                }
                int non_mut_x_u = x_u - xbar_u;
                IntVector choices;
                for (int xbar_v = 0; xbar_v <= x_v; ++xbar_v) {
                    int xxbar_v = x_v - xbar_v;
                    if (xxbar_v > non_mut_x_u) continue;
                    if (xbar_v > xbar_u) continue;

                    choices.push_back(xbar_v);
                }

                boost::random::uniform_int_distribution<> unif(0, choices.size() - 1);
                _xbar[v][mutIdx] = choices[unif(g_rng)];
                _ybar[v][mutIdx] = 0;
            }
        } else {
            if (!(ybar_u > 0)) {
                throw std::runtime_error("Error in samplemutation: ybar_u must be > 0");
            }
            if (y_u == y_v) {
                _xbar[v][mutIdx] = xbar_u;
                _ybar[v][mutIdx] = ybar_u;
            } else if (ybar_u == y_u) {
                _xbar[v][mutIdx] = 0;
                _ybar[v][mutIdx] = y_v;
            } else if (y_u < y_v) {
                // amplification
                int delta = y_v - y_u;
                IntVector choices;
                for (int l = 0; l <= delta; ++l) {
                    choices.push_back(ybar_u + delta);
                }

                boost::random::uniform_int_distribution<> unif(0, choices.size() - 1);
                _xbar[v][mutIdx] = 0;
                _ybar[v][mutIdx] = choices[unif(g_rng)];;
            } else {
                // deletion
                if (!(y_u > y_v)) {
                    throw std::runtime_error("Error in sampleMutations: y_u must be > y_v");
                }
                int non_mut_y_u = y_u - ybar_u;
                IntVector choices;
                for (int ybar_v = 0; ybar_v <= y_v; ++ybar_v) {
                    int yybar_v = y_v - ybar_v;
                    if (yybar_v > non_mut_y_u) continue;
                    if (ybar_v > ybar_u) continue;
                    choices.push_back(ybar_v);
                }

                boost::random::uniform_int_distribution<> unif(0, choices.size() - 1);
                _xbar[v][mutIdx] = 0;
                _ybar[v][mutIdx] = choices[unif(g_rng)];;
            }
        }

        sampleMutation(v, segmentIdx, mutIdx);
    }
}

void Phylogeny::initD(Node u) {
    _D[u].clear();
    for (OutArcIt a(_T, u); a != lemon::INVALID; ++a) {
        Node v = _T.target(a);
        initD(v);
        _D[u].insert(v);
        _D[u].insert(_D[v].begin(), _D[v].end());
    }
}

void Phylogeny::writeProportionFile(std::ostream &out, std::string &outputProportionFilename, int nrSamples) const {
    Digraph::NodeMap<int> nodeToIndex_updated(_T);
    int idx = 0;
    for (NodeIt v(_T); v != lemon::INVALID; ++v) {
        nodeToIndex_updated[v] = idx++;
    }
    if (outputProportionFilename.compare("") != 0) {
        std::ofstream myfile;
        myfile.open(outputProportionFilename);
        if (!myfile) {
            throw std::runtime_error("Error in creating output file for proportions");
        } else {
            myfile << "Node\tProportionofEachSample\n";
            IntNodeMap nodeId(_T);

            for (NodeIt v(_T); v != lemon::INVALID; ++v) {
                myfile << nodeToIndex_updated[v];
                for (int sampleIdx = 0; sampleIdx < nrSamples; sampleIdx++) {
                    myfile << "\t" << _proportions[v][sampleIdx];
                }
                myfile << "\n";
            }


        }
        myfile.close();
    } else {
        throw std::runtime_error("Please specify proportion file name");
    }
}


void Phylogeny::sampleProportions(int nrSamples, double expPurity, double minProportion) {


    const int nrClusters = _clusterToNode.size();

    // make sure to sample mrca and descendants of every cluster
    DoubleVector purityVector(nrSamples);


    for (int sampleIdx = 0; sampleIdx < nrSamples; ++sampleIdx) {
        double samplePurity = expPurity;
        if (expPurity < 1.) {
            double alpha = expPurity * 100;
            double beta = 100 - alpha;
            sftrabbit::beta_distribution<> betaDist(alpha, beta);

            purityVector[sampleIdx] = betaDist(g_rng);
        } else { //added AH
            purityVector[sampleIdx] = 1.;
        }
    }

    int countNodes = 0;
    // initialize proportions
    for (NodeIt v(_T); v != lemon::INVALID; ++v) {
        countNodes++;
        _proportions[v] = DoubleVector(nrSamples, 0.);
        if (v == _root) {
            for (int sampleIdx = 0; sampleIdx < nrSamples; ++sampleIdx) {
                _proportions[v][sampleIdx] = 1. - purityVector[sampleIdx];
            }
        }
    }

    IntMatrix clusterToSample(nrClusters);
    IntVector sampleVector(nrSamples);
    for (int i = 0; i < nrSamples; ++i) {
        sampleVector[i] = i;
    }

    boost::random::uniform_int_distribution<> unif_samples(1, nrSamples);

    IntMatrix sampleToCluster(nrSamples);
    for (int clusterIdx = 0; clusterIdx < nrClusters; ++clusterIdx) {
        std::shuffle(sampleVector.begin(), sampleVector.end(), g_rng);
        int nr_picked_samples = unif_samples(g_rng);

        clusterToSample[clusterIdx] = IntVector(sampleVector.begin(), sampleVector.begin() + nr_picked_samples);
        for (int sample: clusterToSample[clusterIdx]) {
            sampleToCluster[sample].push_back(clusterIdx);
        }
    }

    boost::random::gamma_distribution<> gamma_dist(1, 1);
    for (int sampleIdx = 0; sampleIdx < nrSamples; ++sampleIdx) {
        double minSampleProportion = std::min(minProportion,
                                              purityVector[sampleIdx] / sampleToCluster[sampleIdx].size());

        NodeSet sampledNodes;
        for (int cloneIdx = 0; cloneIdx < sampleToCluster[sampleIdx].size(); ++cloneIdx) {
            if (_clusterD[cloneIdx].size() > 0) {
                boost::random::uniform_int_distribution<> unif_cluster(0, _clusterD[cloneIdx].size() - 1);
                Node v = _clusterD[cloneIdx][unif_cluster(g_rng)];
                sampledNodes.insert(v);
            }
        }
        for (Node v: sampledNodes) {
            int a = _nodeToIndex[v];
        }
        DoubleVector gamma(sampledNodes.size());

        bool ok = false;
        double sum;
        while (!ok) {
            ok = true;
            sum = 0;
            for (int nodeIdx = 0; nodeIdx < sampledNodes.size(); nodeIdx++) {
                gamma[nodeIdx] = gamma_dist(g_rng);
                sum += gamma[nodeIdx];
            }


            for (int nodeIdx = 0; nodeIdx < sampledNodes.size(); nodeIdx++) {
                double prop = gamma[nodeIdx] / sum * purityVector[sampleIdx];
                if (prop < minSampleProportion) ok = false;
            }

        }

        int nodeIdx = 0;
        for (Node n: sampledNodes) {
            double prop = gamma[nodeIdx] / sum * purityVector[sampleIdx];
            _proportions[n][sampleIdx] = prop;
            nodeIdx++;
        }

    }
}


void Phylogeny::initClusterD() {
    for (int clusterIdx = 0; clusterIdx < _clusterToNode.size(); clusterIdx++) {
        initClusterD_helper(_clusterToNode[clusterIdx], clusterIdx);
    }
}

void Phylogeny::initClusterD_helper(Node v, int clusterIdx) {
    IntSet mutations = _clusterToMut[clusterIdx];
    //for (int m : mutations) {
    //    if (_xbar[v][m] == 0 && _ybar[v][m == 0]) {
    //        return;
    //    }
    //}
    _clusterD[clusterIdx].push_back(v);
    for (OutArcIt a(_T, v); a != lemon::INVALID; ++a) {
        initClusterD_helper(_T.target(a), clusterIdx);
    }
}

std::ostream &operator<<(std::ostream &out, const Phylogeny &T) {
    out << T._cnaTrees;
    out << lemon::countNodes(T._T) << " # nodes" << std::endl;

    IntNodeMap nodeId(T._T);
    int index = 0;
    for (NodeIt v(T._T); v != lemon::INVALID; ++v) {
        nodeId[v] = index++;
    }

    out << T._mutToCluster.size() << " #mutations" << std::endl;
    out << T._proportions[T._root].size() << " #samples" << std::endl;

    for (NodeIt v(T._T); v != lemon::INVALID; ++v) {
        out << nodeId[v] << " " << (T._root == v ? "1" : "0")
            << " " << (T._mrca == v ? "1" : "0")
            << " " << (T._trunk[v] ? "1" : "0");

        for (int xbar: T._xbar[v]) {
            out << " " << xbar;
        }

        for (int ybar: T._ybar[v]) {
            out << " " << ybar;
        }

        for (int state: T._charState[v]) {
            out << " " << state;
        }
        for (double prop: T._proportions[v]) {
            out << " " << prop;
        }
        out << std::endl;
    }

    for (ArcIt a(T._T); a != lemon::INVALID; ++a) {
        out << nodeId[T._T.source(a)] << " " << nodeId[T._T.target(a)] << std::endl;
    }

    out << T._clusterToNode.size() << " # mutation clusters" << std::endl;
    for (int clusterIdx = 0; clusterIdx < T._clusterToNode.size(); ++clusterIdx) {
        out << clusterIdx << " " << nodeId[T._clusterToNode[clusterIdx]];
        for (int mutation: T._clusterToMut[clusterIdx]) {
            out << " " << mutation;
        }
        out << std::endl;
    }

    out << T._segmentToMut.size() << " # segments" << std::endl;
    for (int segmentIdx = 0; segmentIdx < T._segmentToMut.size(); ++segmentIdx) {
        out << segmentIdx;
        for (int mutation: T._segmentToMut[segmentIdx]) {
            out << " " << mutation;
        }
        out << std::endl;
    }

    return out;
}

std::istream &operator>>(std::istream &in, Phylogeny &T) {
    in >> T._cnaTrees;

    std::string line;
    getline(in, line);
    std::stringstream ss(line);

    int nrNodes = -1;
    ss >> nrNodes;
    if (nrNodes < 0) {
        throw std::runtime_error("Expected nonnegative number of nodes");
    }

    int nrMutations = -1;
    getline(in, line);
    ss.clear();
    ss.str(line);
    ss >> nrMutations;
    if (nrMutations < 0) {
        throw std::runtime_error("Expected nonnegative number of mutations");
    }

    int nrSamples = -1;
    getline(in, line);
    ss.clear();
    ss.str(line);
    ss >> nrSamples;
    if (nrSamples < 0) {
        throw std::runtime_error("Expected nonnegative number of samples");
    }

    T._T.clear();
    NodeVector indexToNode(nrNodes, lemon::INVALID);
    for (int i = 0; i < nrNodes; ++i) {
        getline(in, line);
        StringVector s;
        boost::split(s, line, boost::is_any_of(" \t"));

        if (s.size() != 1 + 3 + 2 * nrMutations + T._cnaTrees.size() + nrSamples) {
            throw std::runtime_error("Unexpected number of values");
        }

        int nodeIdx = boost::lexical_cast<int>(s[0]);
        Node v = T._T.addNode();
        indexToNode[nodeIdx] = v;
        T._nodeToCluster[v] = -1;
        if (boost::lexical_cast<int>(s[1]) == 1) {
            T._root = v;
        }
        if (boost::lexical_cast<int>(s[2]) == 1) {
            T._mrca = v;
        }
        T._trunk[v] = boost::lexical_cast<int>(s[3]) == 1;

        int offset = 4;
        T._xbar[v] = IntVector(nrMutations, -1);
        for (int mutationIdx = 0; mutationIdx < nrMutations; ++mutationIdx) {
            int xbar = boost::lexical_cast<int>(s[offset + mutationIdx]);
            T._xbar[v][mutationIdx] = xbar;
        }

        offset += nrMutations;
        T._ybar[v] = IntVector(nrMutations, -1);
        for (int mutationIdx = 0; mutationIdx < nrMutations; ++mutationIdx) {
            int ybar = boost::lexical_cast<int>(s[offset + mutationIdx]);
            T._ybar[v][mutationIdx] = ybar;
        }

        offset += nrMutations;

        T._charState[v] = IntVector(T._cnaTrees.size(), -1);
        for (int segmentIdx = 0; segmentIdx < T._cnaTrees.size(); ++segmentIdx) {
            int state = boost::lexical_cast<int>(s[offset + segmentIdx]);
            if (!(0 <= state && state < T._cnaTrees[segmentIdx].k())) {
                throw std::runtime_error("Unexpected state");
            }

            T._charState[v][segmentIdx] = state;
        }

        offset += T._cnaTrees.size();

        T._proportions[v] = DoubleVector(nrSamples, 0);
        for (int sampleIdx = 0; sampleIdx < nrSamples; ++sampleIdx) {
            T._proportions[v][sampleIdx] = boost::lexical_cast<double>(s[offset + sampleIdx]);
        }
    }

    for (int i = 0; i < nrNodes - 1; ++i) {
        int sourceIdx = -1;
        int targetIdx = -1;
        getline(in, line);
        ss.clear();
        ss.str(line);
        ss >> sourceIdx >> targetIdx;

        if (!(0 <= sourceIdx && sourceIdx < nrNodes)) {
            throw std::runtime_error("Invalid edge");
        }

        if (!(0 <= targetIdx && targetIdx < nrNodes)) {
            throw std::runtime_error("Invalid edge");
        }

        T._T.addArc(indexToNode[sourceIdx], indexToNode[targetIdx]);
    }

    int nrClusters = -1;
    getline(in, line);
    ss.clear();
    ss.str(line);
    ss >> nrClusters;

    if (nrClusters < 0) {
        throw std::runtime_error("Invalid number of mutation clusters");
    }

    T._clusterToMut = IntSetVector(nrClusters);
    T._mutToCluster = IntVector(nrMutations);
    T._clusterToNode = NodeVector(nrClusters);
    for (int clusterIdx = 0; clusterIdx < nrClusters; ++clusterIdx) {
        getline(in, line);
        StringVector s;
        boost::split(s, line, boost::is_any_of(" \t"));

        if (s.size() < 2) {
            throw std::runtime_error("Expected at least two values");
        }

        if (boost::lexical_cast<int>(s[0]) != clusterIdx) {
            throw std::runtime_error("Unexpected cluster index");
        }

        int nodeIndex = boost::lexical_cast<int>(s[1]);

        if (!(0 <= nodeIndex && nodeIndex < indexToNode.size())) {
            throw std::runtime_error("Unexpected node index");
        }

        T._clusterToNode[clusterIdx] = indexToNode[nodeIndex];
        T._nodeToCluster[indexToNode[nodeIndex]] = clusterIdx;

        for (int i = 2; i < s.size(); ++i) {
            int mut = boost::lexical_cast<int>(s[i]);
            T._clusterToMut[clusterIdx].insert(mut);
            T._mutToCluster[mut] = clusterIdx;
        }
    }

    int nrSegments = -1;
    getline(in, line);
    ss.clear();
    ss.str(line);
    ss >> nrSegments;

    if (nrSegments != T._cnaTrees.size()) {
        throw std::runtime_error("Invalid number of segments");
    }

    T._segmentToMut = IntSetVector(nrSegments);
    T._mutToSegment = IntVector(nrMutations);

    for (int segmentIdx = 0; segmentIdx < nrSegments; ++segmentIdx) {
        getline(in, line);
        StringVector s;
        boost::split(s, line, boost::is_any_of(" \t"));

        if (s.size() < 1) {
            throw std::runtime_error("Expected at least one value");
        }

        if (boost::lexical_cast<int>(s[0]) != segmentIdx) {
            throw std::runtime_error("Unexpected segment index");
        }

        for (int i = 1; i < s.size(); ++i) {
            int mut = boost::lexical_cast<int>(s[i]);
            T._segmentToMut[segmentIdx].insert(mut);
            T._mutToSegment[mut] = segmentIdx;
        }
    }

    return in;
}

Phylogeny Phylogeny::removeUnsampledNodes() const {
    Phylogeny newPhylo(*this);
    for (ArcIt a(newPhylo._T); a != lemon::INVALID; ++a) {
        Node par = newPhylo._T.source(a);
        int parent = newPhylo._nodeToIndex[par];
        Node ch = newPhylo._T.target(a);
        int child = newPhylo._nodeToIndex[ch];
    }
    for (NodeIt v(newPhylo._T); v != lemon::INVALID; ++v) {
        int b = newPhylo._nodeToIndex[v];
    }

    BoolNodeMap sampled(newPhylo._T, false);
    NodeSet unsampledNodes;
    for (NodeIt v(newPhylo._T); v != lemon::INVALID; ++v) {
        bool sampled_v = false;
        for (double prop: newPhylo._proportions[v]) {
            sampled_v |= prop > 0.;
        }
        sampled[v] = sampled_v;
        if (!sampled_v) {
            unsampledNodes.insert(v);
        }
        if (sampled_v) {
            int b = 5;
        }
    }

    while (true) {
        Node toRemove = lemon::INVALID;
        int outDeg = -1;
        for (Node v: unsampledNodes) {
            outDeg = lemon::countOutArcs(newPhylo._T, v);
            if (outDeg <= 1) {
                toRemove = v;
                break;
            }
        }

        if (toRemove == lemon::INVALID) break;
        //assert(toRemove != newPhylo._root); //AH 10/11
        Node child = lemon::INVALID;
        if (toRemove == newPhylo._root) {
            if (outDeg != 1) {
                throw std::runtime_error("Error: out degree of root being removed is not 1");
            }
            child = newPhylo._T.target(OutArcIt(newPhylo._T, toRemove));
            newPhylo._root = child;
            newPhylo._trunkLength--;
        } else {
            Node parent = newPhylo._T.source(InArcIt(newPhylo._T, toRemove));
            if (outDeg == 1) {
                child = newPhylo._T.target(OutArcIt(newPhylo._T, toRemove));
                newPhylo._T.addArc(parent, child);
            }
            if (newPhylo._trunk[toRemove]) {
                newPhylo._trunkLength--;
            }
            if (toRemove == newPhylo._mrca) {
                if (outDeg != 1) {
                    throw std::runtime_error("Error: out degree of mrca being removed is not 1");
                }
                newPhylo._mrca = child;
            }
        }
        //updating nodeToCluster and clusterToNode assignments for newPhylo
        int assignedCluster = newPhylo._nodeToCluster[toRemove];
        if (assignedCluster != -1) {
            newPhylo._nodeToCluster[child] = assignedCluster;
            newPhylo._clusterToNode[assignedCluster] = child;
        }
        unsampledNodes.erase(toRemove);
        newPhylo._T.erase(toRemove);
    }

    newPhylo.initD(newPhylo._root);
    newPhylo._clusterD = NodeMatrix(newPhylo._clusterToNode.size(), NodeVector(0));
    //newPhylo.initClusterD(newPhylo._mrca, 0);
    newPhylo.initClusterD();
    for (ArcIt a(newPhylo._T); a != lemon::INVALID; ++a) {
        Node par = newPhylo._T.source(a);
        int parent = newPhylo._nodeToIndex[par];
        Node ch = newPhylo._T.target(a);
        int child = newPhylo._nodeToIndex[ch];
    }


    return newPhylo;
}