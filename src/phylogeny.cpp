//
// Created by Mohammed El-Kebir on 8/20/23.
//

#include "phylogeny.h"
#include <algorithm>
#include <lemon/connectivity.h>

Phylogeny::Phylogeny()
  : _T()
  , _root(_T.addNode())
  , _mrca(lemon::INVALID)
  , _trunk(_T)
  , _cnaTrees()
  , _charState(_T)
  , _trunkLength(0)
  , _mutToCluster()
  , _clusterToMut()
  , _xbar(_T)
  , _ybar(_T)
  , _clusterToNode()
  , _nodeToCluster(_T)
  , _mutToSegment()
  , _segmentToMut()
  , _D(_T)
{
  _trunk[_root] = false;
}

void Phylogeny::addSegment(const CnaTree& cnaTree, bool truncal)
{
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

void Phylogeny::update(const CnaTree& cnaTree, bool truncal,
                       const Digraph& G, const NodePairNodeMap& mapG, const NodeNodeSetMap& toG,
                       const NodeNodeSetMap& cnaToG, const Node rootG, const SubDigraph& subG)
{
  NodeSet toErase;
  for (NodeIt vv(_T); vv != lemon::INVALID; ++vv)
  {
    toErase.insert(vv);
  }

  NodeNodeMap productNodesToNewNodes(G, lemon::INVALID);
  for (SubNodeIt uv(subG); uv != lemon::INVALID; ++uv)
  {
    if (SubOutArcIt(subG, uv) != lemon::INVALID || SubInArcIt(subG, uv) != lemon::INVALID)
    {
      Node new_uv = _T.addNode();
      _nodeToCluster[new_uv] = -1;
      _charState[new_uv] = _charState[mapG[uv].first];
      _charState[new_uv].push_back(cnaTree.state(mapG[uv].second));
      _trunk[new_uv] = false;
//      _trunk[new_uv] = _trunk[mapG[uv].first];
      productNodesToNewNodes[uv] =  new_uv;

      if (mapG[uv].first == _root && cnaTree.state(mapG[uv].second) == cnaTree.rootState())
      {
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

  for (SubArcIt aa(subG); aa != lemon::INVALID; ++aa)
  {
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
  while (lemon::countOutArcs(_T, _mrca) == 1)
  {
    _trunk[_mrca] = true;
    _mrca = _T.target(OutArcIt(_T, _mrca));
    ++trunkLength;
    if (trunkLength == _trunkLength + 1)
      break;
  }
  _trunkLength = trunkLength;

  for (Node node : toErase)
  {
    _T.erase(node);
  }
}

void Phylogeny::sample(const Digraph& G, const NodePairNodeMap& mapG, const NodeNodeSetMap& toG,
                       const NodeNodeSetMap& cnaToG, const Node rootG, SubDigraph& subG) const
{
  // build frontier
  ArcSet F;
  for (OutArcIt a(G, rootG); a != lemon::INVALID; ++a)
  {
    F.insert(a);
  }

  ArcSet sampledEdges;
  while (!F.empty())
  {
    // pick edge at random
    Arc a = lemon::INVALID;

    int sample = std::uniform_int_distribution<>(1, F.size())(g_rng);
    int i = 1;
    for (Arc aa : F)
    {
      if (i++ == sample)
      {
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

    if (mapG[uv_source].first == mapG[uv_target].first)
    {
      // remove all incoming arcs to mapG[uv_target].second
      for (Node uv : cnaToG[mapG[uv_target].second])
      {
        for (InArcIt in_uv(G, uv); in_uv != lemon::INVALID; ++in_uv)
        {
          if (sampledEdges.count(in_uv) == 0 && (mapG[G.source(in_uv)].second == mapG[uv_source].second) && (mapG[G.target(in_uv)].second == mapG[uv_target].second))
          {
            subG.disable(in_uv);
            F.erase(in_uv);
          }
        }
      }
    }

    if (mapG[uv_source].second == mapG[uv_target].second)
    {
      // remove all incoming arcs to mapG[uv_target].first
      for (Node uv : toG[mapG[uv_target].first])
      {
        for (InArcIt in_uv(G, uv); in_uv != lemon::INVALID; ++in_uv)
        {
          if (sampledEdges.count(in_uv) == 0 && (mapG[G.source(in_uv)].first == mapG[uv_source].first) && (mapG[G.target(in_uv)].first == mapG[uv_target].first))
          {
            subG.disable(in_uv);
            F.erase(in_uv);
          }
        }
      }
    }

    for (SubOutArcIt aa(subG, uv_target); aa != lemon::INVALID; ++aa)
    {
      F.insert(aa);
    }
  }
}

void Phylogeny::createProductGraph(const CnaTree& cnaTree, Digraph& G, NodePairNodeMap& mapG,
                                   NodeNodeSetMap& toG, NodeNodeSetMap& cnaToG, Node& rootG) const
{
  // add nodes to product graph
  for (NodeIt u(_T); u != lemon::INVALID; ++u)
  {
    for (NodeIt v(cnaTree.T()); v != lemon::INVALID; ++v)
    {
      Node uv = G.addNode();
      mapG[uv] = std::make_pair(u, v);
      toG[u].insert(uv);
      cnaToG[v].insert(uv);

      if (u == _root && cnaTree.state(v) == cnaTree.rootState())
      {
        rootG = uv;
      }
    }
  }

  // add edges to product graph
  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    for (Node uv : toG[_T.source(a)])
    {
      for (Node uv2 : toG[_T.target(a)])
      {
        if (mapG[uv].second == mapG[uv2].second)
        {
          G.addArc(uv, uv2);
        }
      }
    }
  }
  for (ArcIt a(cnaTree.T()); a != lemon::INVALID; ++a)
  {
    for (Node uv : cnaToG[cnaTree.T().source(a)])
    {
      for (Node uv2 : cnaToG[cnaTree.T().target(a)])
      {
        if (mapG[uv].first == mapG[uv2].first && !_trunk[mapG[uv].first])
        {
          G.addArc(uv, uv2);
        }
      }
    }
  }
}

void Phylogeny::printProductGraph(const Digraph& G,
                                  const SubDigraph& subG, const NodePairNodeMap& mapG, const CnaTree& cnaTree,
                                  std::ostream& out) const
{
  out << "digraph G {" << std::endl;

  for (SubNodeIt v(subG); v != lemon::INVALID; ++v)
  {
    out << "\t" << subG.id(v) << " [label=\"";
    int idx = 0;
    for (int state : _charState[mapG[v].first])
    {
      out << _cnaTrees[idx].label(state) << ",";
      idx++;
    }
    out << cnaTree.label(cnaTree.state(mapG[v].second));
    if (mapG[v].first == _root)
    {
      out << "\\nroot";
    }
    if (mapG[v].first == _mrca)
    {
      out << "\\nmrca";
    }
    if (_trunk[mapG[v].first])
    {
      out << "\\ntrunk";
    }
    out << "\"";
    if (!(SubOutArcIt(subG, v) != lemon::INVALID || SubInArcIt(subG, v) != lemon::INVALID))
    {
      out << ",style=dashed";
    }
    if (_trunk[mapG[v].first])
    {
      out << ",color=red";
    }
    out << "]" << std::endl;
  }

  for (ArcIt a(G); a != lemon::INVALID; ++a)
  {
    out << "\t" << G.id(G.source(a)) << " -> " << G.id(G.target(a));
    if (!subG.status(a))
    {
      out << " [style=dashed]";
    }
    out << std::endl;
  }

  out << "}" << std::endl;
}

void Phylogeny::writeDOT(std::ostream& out) const
{
  out << "digraph T {"<< std::endl;
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    out << "\t" << _T.id(v) << " [label=\"";
    for (int i = 0; i < _charState[v].size(); ++i)
    {
//      out << _cnaTrees[i].label(_charState[v][i]) << " ; ";
      out << _cnaTrees[i].label(_charState[v][i]) << " ; ";
    }
    for (int segmentIdx = 0; segmentIdx < _cnaTrees.size(); ++segmentIdx)
    {
      out << "\\n";
      for (int mutIdx : _segmentToMut[segmentIdx])
      {
        out << _xbar[v][mutIdx] << "|" << _ybar[v][mutIdx] << " ; ";
      }
    }
    if (v == _root)
    {
      out << "\\nroot";
    }
    if (v == _mrca)
    {
      out << "\\nmrca";
    }
    if (_trunk[v])
    {
      out << "\\ntrunk";
    }
    if (_nodeToCluster[v] != -1)
    {
      out << "\\ncluster " << _nodeToCluster[v];
    }
    out << "\"]" << std::endl;
  }

  for (ArcIt a(_T); a != lemon::INVALID; ++a)
  {
    out << _T.id(_T.source(a)) << " -> " << _T.id(_T.target(a)) << std::endl;
  }

  out << "}" << std::endl;
}

void Phylogeny::sampleMutations(int n, int l)
{
  assert(n >= 1);
  assert(l >= 1);
  assert(l <= n);

  const int k = getNrSegments();
  _clusterToMut = IntSetVector(l);

  // assign mutations to segments
  std::uniform_int_distribution<> uniform_segments(0, k - 1);

  for (int i = 0; i < n; ++i)
  {
    int segmentIdx = uniform_segments(g_rng);
    _mutToSegment.push_back(segmentIdx);
    _segmentToMut[segmentIdx].insert(i);
  }

  // assign mutations
  std::uniform_int_distribution<> uniform_clusters(0, l - 1);

  for (int i = 0; i < n; ++i)
  {
    int clusterIdx = uniform_clusters(g_rng);
    _mutToCluster.push_back(clusterIdx);
    _clusterToMut[clusterIdx].insert(i);
  }

  NodeVector trunkNodes;
  trunkNodes.push_back(_mrca);
  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    if (_trunk[v] && v != _root)
    {
      trunkNodes.push_back(v);
    }
  }

  if (trunkNodes.empty())
  {
    trunkNodes.push_back(_root);
  }

  // pick trunk mutation cluster
  std::uniform_int_distribution<> unif1(0, trunkNodes.size() - 1);
  int idx = unif1(g_rng);
  _clusterToNode.push_back(trunkNodes[idx]);
  _nodeToCluster[trunkNodes[idx]] = 0;

  // sample the remaining l-1 mutation cluster locations
  NodeVector remainingNodes(_D[_clusterToNode[0]].begin(), _D[_clusterToNode[0]].end());
  std::shuffle(remainingNodes.begin(), remainingNodes.end(), g_rng);
  assert(remainingNodes.size() >= l - 1);

  for (int i = 0; i < l - 1; ++i)
  {
    _clusterToNode.push_back(remainingNodes[i]);
    _nodeToCluster[remainingNodes[i]] = i + 1;
  }

  std::uniform_int_distribution<> unif01(0,1);

  for (NodeIt v(_T); v != lemon::INVALID; ++v)
  {
    _xbar[v] = IntVector(n, 0);
    _ybar[v] = IntVector(n, 0);
  }

  // sample xbar and ybar for each mutation
  for (int mutIdx = 0; mutIdx < n; ++mutIdx)
  {
    int clusterIdx = _mutToCluster[mutIdx];
    int segmentIdx = _mutToSegment[mutIdx];
    Node mutationNode = _clusterToNode[clusterIdx];

    int x = _cnaTrees[segmentIdx].x(_charState[mutationNode][segmentIdx]);
    int y = _cnaTrees[segmentIdx].y(_charState[mutationNode][segmentIdx]);

    if (x > 0 && y > 0)
    {
      bool mut_x = unif01(g_rng) == 1;
      if (mut_x)
      {
        _xbar[mutationNode][mutIdx] = 1;
        _ybar[mutationNode][mutIdx] = 0;
        sampleMutation(mutationNode, segmentIdx, mutIdx);
      }
      else
      {
        _xbar[mutationNode][mutIdx] = 0;
        _ybar[mutationNode][mutIdx] = 1;
        sampleMutation(mutationNode, segmentIdx, mutIdx);
      }
    }
    else if (x > 0)
    {
      _xbar[mutationNode][mutIdx] = 1;
      _ybar[mutationNode][mutIdx] = 0;
      sampleMutation(mutationNode, segmentIdx, mutIdx);
    }
    else if (y > 0)
    {
      _xbar[mutationNode][mutIdx] = 0;
      _ybar[mutationNode][mutIdx] = 1;
      sampleMutation(mutationNode, segmentIdx, mutIdx);
    }
  }
}

void Phylogeny::sampleMutation(const Node u, const int segmentIdx, const int mutIdx)
{
  int x_u = _cnaTrees[segmentIdx].x(_charState[u][segmentIdx]);
  int y_u = _cnaTrees[segmentIdx].y(_charState[u][segmentIdx]);
  int xbar_u = _xbar[u][mutIdx];
  int ybar_u = _ybar[u][mutIdx];
  assert(0 <= xbar_u && xbar_u <= x_u);
  assert(0 <= ybar_u && ybar_u <= y_u);

  for (OutArcIt a(_T, u); a != lemon::INVALID; ++a)
  {
    Node v = _T.target(a);
    int x_v = _cnaTrees[segmentIdx].x(_charState[v][segmentIdx]);
    int y_v = _cnaTrees[segmentIdx].y(_charState[v][segmentIdx]);

    if (xbar_u == 0 && ybar_u == 0)
    {
      _xbar[v][mutIdx] = _ybar[v][mutIdx];
    }
    else if (xbar_u > 0)
    {
      if (x_u == x_v)
      {
        _xbar[v][mutIdx] = xbar_u;
        _ybar[v][mutIdx] = ybar_u;
      }
      else if (xbar_u == x_u)
      {
        _xbar[v][mutIdx] = x_v;
        _ybar[v][mutIdx] = 0;
      }
      else if (x_u < x_v)
      {
        // amplification
        int delta = x_v - x_u;
        IntVector choices;
        for (int l = 0; l <= delta; ++l)
        {
          choices.push_back(xbar_u + delta);
        }

        std::uniform_int_distribution<> unif(0, choices.size() - 1);
        _xbar[v][mutIdx] = choices[unif(g_rng)];
        _ybar[v][mutIdx] = 0;
      }
      else
      {
        // deletion
        assert(x_u > x_v);
        int non_mut_x_u = x_u - xbar_u;
        IntVector choices;
        for (int xbar_v = 0; xbar_v <= x_v; ++xbar_v)
        {
          int xxbar_v = x_v - xbar_v;
          if (xxbar_v > non_mut_x_u) continue;
          if (xbar_v > xbar_u) continue;

          choices.push_back(xbar_v);
        }

        std::uniform_int_distribution<> unif(0, choices.size() - 1);
        _xbar[v][mutIdx] = choices[unif(g_rng)];
        _ybar[v][mutIdx] = 0;
      }
    }
    else
    {
      assert(ybar_u > 0);
      if (y_u == y_v)
      {
        _xbar[v][mutIdx] = xbar_u;
        _ybar[v][mutIdx] = ybar_u;
      }
      else if (ybar_u == y_u)
      {
        _xbar[v][mutIdx] = 0;
        _ybar[v][mutIdx] = y_v;
      }
      else if (y_u < y_v)
      {
        // amplification
        int delta = y_v - y_u;
        IntVector choices;
        for (int l = 0; l <= delta; ++l)
        {
          choices.push_back(ybar_u + delta);
        }

        std::uniform_int_distribution<> unif(0, choices.size() - 1);
        _xbar[v][mutIdx] = 0;
        _ybar[v][mutIdx] = choices[unif(g_rng)];;
      }
      else
      {
        // deletion
        assert(y_u > y_v);
        int non_mut_y_u = y_u - ybar_u;
        IntVector choices;
        for (int ybar_v = 0; ybar_v <= y_v; ++ybar_v)
        {
          int yybar_v = y_v - ybar_v;
          if (yybar_v > non_mut_y_u) continue;
          if (ybar_v > ybar_u) continue;
          choices.push_back(ybar_v);
        }

        std::uniform_int_distribution<> unif(0, choices.size() - 1);
        _xbar[v][mutIdx] = 0;
        _ybar[v][mutIdx] = choices[unif(g_rng)];;
      }
    }

    sampleMutation(v, segmentIdx, mutIdx);
  }
}

void Phylogeny::initD(Node u)
{
  for (OutArcIt a(_T, u); a != lemon::INVALID; ++a)
  {
    Node v = _T.target(a);
    initD(v);
    _D[u].insert(v);
    _D[u].insert(_D[v].begin(), _D[v].end());
  }
}