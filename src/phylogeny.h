//
// Created by Mohammed El-Kebir on 8/20/23.
//

#ifndef PHYLOGENY_H
#define PHYLOGENY_H

#include "utils.h"
#include "cnatree.h"

class Phylogeny
{
public:
  /// Default constructor
  Phylogeny();

  /// Copy constructor
  Phylogeny(const Phylogeny& other);

  /// Add segment
  void addSegment(const CnaTree& cnaTree, bool truncal);

  /// Sample mutations
  ///
  /// \param n Number of mutations
  /// \param l Number of mutation clusters
  void sampleMutations(int n, int l);

  /// Write phylogeny in graphviz format
  ///
  /// \param out Output stream
  void writeDOT(std::ostream& out) const;

  /// Sample mixture proportions on nodes
  ///
  /// \param nrSamples Number of samples
  /// \param expPurity Expected purity
  /// \param minProportion Minimum proportion
  void sampleProportions(int nrSamples, double expPurity, double minProportion);

  int getNrSegments() const
  {
    return _cnaTrees.size();
  };

  int getTrunkLength() const
  {
    return _trunkLength;
  };

  Node getRoot() const
  {
    return _root;
  }

  Node getMRCA() const
  {
    return _mrca;
  }

  Phylogeny removeUnsampledNodes() const;

private:
  typedef Digraph::NodeMap<IntVector> CharacterStateVector;
  typedef Digraph::NodeMap<std::pair<Node, Node> > NodePairNodeMap;

  void createProductGraph(const CnaTree& cnaTree, Digraph& G, NodePairNodeMap& mapG,
                          NodeNodeSetMap& toG, NodeNodeSetMap& cnaToG, Node& rootG) const;

  void sample(const Digraph& G, const NodePairNodeMap& mapG,
              const NodeNodeSetMap& toG, const NodeNodeSetMap& cnaToG, const Node rootG,
              SubDigraph& subG) const;

  void update(const CnaTree& cnaTree, bool truncal,
              const Digraph& G, const NodePairNodeMap& mapG, const NodeNodeSetMap& toG,
              const NodeNodeSetMap& cnaToG, const Node rootG, const SubDigraph& subG);

  void printProductGraph(const Digraph& G,
                         const SubDigraph& subG, const NodePairNodeMap& mapG, const CnaTree& cnaTree,
                         std::ostream& out) const;

  void initD(Node v);

  void initClusterD(Node v, int clusterIdx);

  void sampleMutation(const Node u, const int segmentIdx, const int mutIdx);

private:
  /// Tree
  Digraph  _T;
  /// Root
  Node _root;
  /// MRCA
  Node _mrca;
  /// Truncal indicator, if trunk[v] == true then v has at most 1 child
  BoolNodeMap _trunk;
  /// Cna tree vector
  CnaTreeVector _cnaTrees;
  /// Character state vector
  CharacterStateVector _charState;
  /// Trunk length
  int _trunkLength;
  /// Mutation to cluster
  IntVector _mutToCluster;
  /// Cluster to mutations
  IntSetVector _clusterToMut;
  /// Number of mutated maternal copies
  IntVectorNodeMap _xbar;
  /// Number of mutated paternal copies
  IntVectorNodeMap _ybar;
  /// Cluster to node assignment
  NodeVector _clusterToNode;
  /// Node to cluster assignment
  IntNodeMap _nodeToCluster;
  /// Mutation to segment assignment
  IntVector _mutToSegment;
  /// Segment to mutation assignment
  IntSetVector _segmentToMut;
  /// Descendant set
  NodeNodeSetMap _D;
  /// Cluster descendants
  NodeMatrix _clusterD;
  /// Proportions
  DoubleVectorNodeMap _proportions;

  friend std::ostream& operator<<(std::ostream& out, const Phylogeny& T);
  friend std::istream& operator>>(std::istream& in, Phylogeny& T);
};

#endif //PHYLOGENY_H
