//
// Created by Leah Weber on 5/16/24.
//

#ifndef CLONALTREE_H
#define CLONALTREE_H

#include "utils.h"
#include "basetree.h"
#include "genotypetree.h"
#include "cnatree.h"
#include <memory>

struct Data
{
    /// Default constructor
    Data()
            : _Var()
            , _Total()
            , _X()
            , _Y()
    {
    }

    /// Constructor
    Data(IntMatrix& var, IntMatrix& total, IntMatrix& xobs, IntMatrix& yobs)
            : _Var(var)
            , _Total(total)
            , _X(xobs)
            , _Y(yobs)
    {
    }


    /// Variant read count matrix n x m (cells x snvs)
    IntMatrix _Var;
    /// Total read count matrix n x m (cells x snvs)
    IntMatrix _Total;
    /// Number of  maternal copies matrix (cells x segments)
    IntMatrix _X;
    /// Number of  paternal copies matrix (cells x segments)
    IntMatrix _Y;
};

typedef GenotypeTree::Genotype Genotype;
typedef CnaTree::CnaGenotype CnaGenotype;
typedef std::map<int, Genotype> IntToGenotypeMap;
typedef std::map<int, CnaGenotype> IntToCnaGenotypeMap;
typedef Digraph :: NodeMap<IntToGenotypeMap> Genotypes;
typedef Digraph :: NodeMap<IntToCnaGenotypeMap> CnStates;
typedef std::map<int,int> IntMap;
typedef std::map<int,IntSet> IntSetMap;
typedef std::set<CnaGenotype> CnaGenotypeSet;




class ClonalTree: public BaseTree {

public:
    ClonalTree(); ///default constructor

    typedef Digraph::NodeMap<std::map<int, GenotypeTreeVector>> GenotypeTreeNodeMap;
//    ClonalTree() = default; // Default constructor
//
//     Delete the copy constructor and copy assignment operator
//    ClonalTree(const ClonalTree&) = delete;
//    ClonalTree& operator=(const ClonalTree&) = delete;
////
////     Provide move constructor and move assignment operator
//    ClonalTree(ClonalTree&&) noexcept = default;
//    ClonalTree& operator=(ClonalTree&&) noexcept = default;

//    ClonalTree(const ClonalTree& other);

    ClonalTree(const std::vector<IntPair>& edges,
               const std::map<int, IntToGenotypeMap>& genotypes,
               const std::map<int, IntToCnaGenotypeMap>& cnStates,
               const IntMap mut2seg, const IntSet& mutClusts);

    const IntMap& getNodes() const {
        return _nodes;
    }

    double computeLikelihood(const IntSetMap& phi,
                            const Data& D,
                             double alpha,
                             double lambda);



    CnaTree getCNATree(int segment);
    GenotypeTree getGenotypeTree(int snv);

    /// given fixed genotypes, a clonal tree optimize the cell assignments given the Data
    std::pair<double, IntSetMap>  assignCells(const IntSet& cells,
                                                          const Data& D,
                                                          double alpha,
                                                          double lambda);

    /// given a cell mapping phi, optimize the genotypes of SNVs given the Tree
    void assignGenotypes(const IntSetMap& phi, const Data& D,double alpha, double lambda);

    /// jointly find the optimal cell mapping phi and genotypes given the Data and clonal tree T
    std::pair<double, IntSetMap> optimize(const IntSet &cells, const Data &D, double alpha, double lambda);

//    void genotypeTreeToClonalTree( int mutClustIdx,  const GenotypeTree& T, const IntSet& muts, int seg);
    void computeSnvAssignmentCost(int mutClustIdx, int treeIdx, std::map<int, double> &snvCost,
                                              std::map<int, std::pair<int, int>> &snvAssign,
                                              const GenotypeTree& T, const IntSetMap& phi, const IntSet& muts,
                                              int seg, const Data& data, double alpha, double lambda);

    IntSetMap collapsePhi(const ClonalTree& ct, const GenotypeTree& T, const IntSetMap& phi, const IntSet& descNodeStates, int seg);
private:
    Genotypes _genotypes;
    CnStates _cnStates;
    IntMap _mut2seg;
    IntSetMap _seg2muts;
    std::map<int,int> _nodes;
    IntSet _segments;
    IntSet _muts;
    BoolNodeMap _isMutCluster;

    GenotypeTreeNodeMap  _rho;

    double nodeCNAcost(const Node& u, const IntSet& cells, const IntMatrix& xobs, const IntMatrix& yobs);
    double nodeSNVcost(const Node& u, const IntSet& cells,const IntMatrix& var, const IntMatrix& total, double alpha);
    static double binomialProbability(int n, int k, double p);
    CnaGenotypeSet descCNAGenotypes(const Node& u, int segment);

    void updateGenotype(int mut, Node u, GenotypeTree T);
    void updateGenotype(IntSet muts, Node u, GenotypeTree T);

    double computeGenotypeLikelihood(int j,IntSetMap phi, const Data& D, double alpha, double lambda);


};


#endif //CLONALTREE_H
