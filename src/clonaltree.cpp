//
// Created by Leah Weber on 5/16/24.
//

#include "clonaltree.h"
#include "utils.h"
#include <lemon/connectivity.h>
#include <boost/math/distributions/binomial.hpp>

ClonalTree::ClonalTree(const std::vector<IntPair>& edges,
                       const std::map<int, IntToGenotypeMap>& genotypes,
                       const std::map<int, IntToCnaGenotypeMap>& cnStates,
                       const IntMap mut2seg, const IntSet& mutClusts)
        : BaseTree(),
          _genotypes(_T),
          _cnStates(_T),
          _isMutCluster(_T, false),
          _rho(_T),
          _mut2seg(mut2seg),
          _seg2muts(),
          _segments(),
          _nodes()
{
    //create the mapping of segments to SNVs
    for(const auto& pair : _mut2seg)
    {
        int mut = pair.first;
        int seg = pair.second;
        _seg2muts[seg].insert(mut);
        _segments.insert(seg);
        _muts.insert(mut);
    }
    //initialize the BaseTree

    _k = edges.size() + 1;
    for(const auto& edge: edges)
    {
        _nodes[edge.second] = _nodes.size() + 1;
    }

    //find the root
    int root;
    for (const auto& edge: edges)
    {
        if (_nodes.count(edge.first) == 0)
        {
            root = edge.first;
            break;
        }
    }
    _nodes[root] = 0;
    _stateToNode = NodeVector(_nodes.size(), lemon::INVALID);

    // 3. add nodes
    char buf[1024];
    for (const auto& node: _nodes)
    {
        Node v = _T.addNode();
        _nodeToState[v] =node.second;
        _stateToNode[node.second] = v;

        if(mutClusts.count(node.first) > 0)
        {
            _isMutCluster[v] = true;
        }


        // Assign genotypes and cnStates
        if (genotypes.find(node.first) != genotypes.end()) {
            _genotypes[v] = genotypes.at(node.first);
        }

        if (cnStates.find(node.first) != cnStates.end()) {
            _cnStates[v] = cnStates.at(node.first);
        }
        snprintf(buf, 1024, "%d", node.first);
        _label[v] = buf;
    }

    _root = _stateToNode[0];

    // 4. add edges
    for (const auto& edge: edges)
    {
        _T.addArc(_stateToNode[_nodes[edge.first]], _stateToNode[_nodes[edge.second]]);
    }



    initD(_root);

    //5. compute rho, a mapping of SNV  cluster to valid SNV trees
    for(int ell: _segments)
    {
        CnaTree S = getCNATree(ell);
//        std::cout << S << std::endl;


            for(NodeIt u(_T); u != lemon::INVALID; ++u )
            {

                if(_isMutCluster[u]){
//                    std::cout<< _nodeToState[u] << std::endl;
                    CnaGenotype uCnaGeno = _cnStates[u][ell];

                    CnaGenotypeSet descCNAs =descCNAGenotypes(u,ell);
                    GenotypeTree::GenotypeEdgeSetSet snvTrees;
                    S.enumerateGenotypeTrees(snvTrees);
                    for(auto edgeSet: snvTrees)
                    {

                        GenotypeTree T(edgeSet);
//                        std::cout << T << std::endl;
                        int mutState = T.mutationState();
                        ///check that the split state of SNV Tree T matches the copy number state of Node u
                        if(T.x(mutState) == uCnaGeno._x && T.y(mutState) == uCnaGeno._y)
                        {
                            IntSet descMutState =T.D(mutState);
                            CnaGenotypeSet mutTreeCnaDesc;
                            for(int j: descMutState)
                            {
                                CnaGenotype g(T.x(j), T.y(j));
                                if(g != uCnaGeno){
                                    mutTreeCnaDesc.insert(g);
                                }

                            }

                            if(mutTreeCnaDesc == descCNAs)
                            {
                                _rho[u][ell].push_back(T);
                            }

                        }

                    }


            }


        }
    }
    assert(lemon::dag(_T));


}

CnaGenotypeSet ClonalTree::descCNAGenotypes(const Node& u,int ell )
{
    CnaGenotypeSet descCNAGenos;
    CnaGenotype genoU = _cnStates[u][ell];
    IntSet desc = D(_nodeToState[u]);
    double cost = 0.0;

    for(int i: desc){
       if(genoU != _cnStates[node(i)][ell])
       {
           descCNAGenos.insert(_cnStates[node(i)][ell]);
       }
    }
   return descCNAGenos;
}

double ClonalTree::assignGenotypes(const IntSetMap& phi, const Data& data,
                                 double alpha, double lambda)
{
    ///for every snv j, try to assign it to each snv cluster q and snv tree _rho[q]
    double cost = 0.0;
    for(int ell: _segments)
    {
        ///tracking snv to cluster and tree assignment;
        std::map<int, std::pair<int, int>> snvAssign;
        std::map<int,double> snvCost;

        ///initialize the SNV likelihood
        for(int j: _seg2muts[ell])
        {
          snvCost[j] = INFINITY;

        }
        for(NodeIt u(_T); u != lemon::INVALID; ++u)
        {
            if(_isMutCluster[u])
            {
                int id = state(u);
                int treeIdx =0;
                for(GenotypeTree T: _rho[u][ell]){
                    computeSnvAssignmentCost(id, treeIdx, snvCost, snvAssign,
                                             T, phi, _seg2muts[ell], ell, data, alpha, lambda);

                    treeIdx++;
                }

            }
            ///Find the best assignment of SNV to SNV cluster and SNV tree and update genotype

        }

        std::map<std::pair<int, int>, IntSet> clusterGenotypeTreeToMutset;
        for(int j: _seg2muts[ell])
        {
            std::pair<int,int> ass =snvAssign[j];
            cost += snvCost[j];
            clusterGenotypeTreeToMutset[ass].insert(j);
            /// to update a single SNV at a time
//            Node n =_stateToNode[snvAssign[j].first];
//            GenotypeTree T = _rho[n][ell][snvAssign[j].second];
//            updateGenotype(j,n,T);
        }

        /// do batch updating of SNVs that are mapped to the same cluster and same SNV tree;
        for(auto clust2tree: clusterGenotypeTreeToMutset)
        {
            int clust = clust2tree.first.first;
            int treeIdx = clust2tree.first.second;
            IntSet muts = clust2tree.second;

            GenotypeTree T = _rho[node(clust)][ell][treeIdx];
            updateGenotype(muts, node(clust), T);

        }





    }

    return cost;
}

void ClonalTree::computeSnvAssignmentCost(int mutClustIdx, int treeIdx, std::map<int, double> &snvCost,
                                          std::map<int, std::pair<int, int>> &snvAssign,
                                          const GenotypeTree& T, const IntSetMap& phi, const IntSet& muts,
                                          int seg, const Data& data, double alpha, double lambda)
{
    std::vector<IntPair> edges;
    std::map<int, IntToGenotypeMap> genotypes;
    std::map<int, IntToCnaGenotypeMap> cnStates;
    IntMap mut2seg;
    IntSet mutClusts;

    CnaGenotype cn = _cnStates[node(mutClustIdx)][seg];
    IntSet descNodeStates = D(mutClustIdx);


    std::cout << T << std::endl;
    int mutState = T.mutationState();

    for(int i=0; i < T.k(); ++i)
    {
        int x;
        int y;
        int xbar;
        int ybar;



        if(i == mutState){
            mutClusts.insert(i);
        }
        T.getCnState(i, x, y,xbar, ybar);

        for(int j: muts)
        {
            genotypes[i][j]=Genotype(x, y, xbar, ybar);

            mut2seg[j] = seg;

        }
        cnStates[i][seg] = CnaGenotype(x, y);


    }


    for(ArcIt a(T.T()); a != lemon::INVALID; ++a)
    {
        int u  = T.state(T.T().source(a));
        int v =  T.state(T.T().target(a));
        edges.push_back(std::make_pair(u,v));
    }

    ClonalTree ct(edges, genotypes, cnStates, mut2seg, mutClusts);
    IntSetMap genoPhi = collapsePhi(ct, T, phi, descNodeStates, seg);


    for(int j: muts)
    {
       double snvLikelihood = ct.computeGenotypeLikelihood(j,genoPhi, data , alpha, lambda);
       if(snvLikelihood <= snvCost[j])
       {
           snvCost[j] = snvLikelihood;
           snvAssign[j] = {mutClustIdx, treeIdx};
       }


    }

}

IntSetMap ClonalTree::collapsePhi(const ClonalTree& ct, const GenotypeTree& T, const IntSetMap& phi, const IntSet& descNodeStates, int seg)
{
    IntSetMap newPhi;
    int newNode;
    CnaGenotype splitCN(T.x(T.mutationState()),T.y(T.mutationState() ));
    std::cout << splitCN._x << "|" << splitCN._y << std::endl;
    for (auto cellpair: phi)
    {

        CnaGenotype cellCn = _cnStates[node(cellpair.first)][seg];
        bool mutPres = descNodeStates.count(cellpair.first);
        newNode = -1;

        for (auto pair: ct.getNodes())
        {
            int x = T.x(pair.first);
            int y = T.y(pair.first);

            CnaGenotype newCn(x, y);

            if (splitCN == newCn && cellCn == newCn)
            {


                if(pair.first == T.mutationState() && mutPres)
                {
                    /// if it's the mutation state there could be two possible causes
                    /// case 1
                   newNode = pair.second;
                   break;



                }else if(!mutPres && pair.first != T.mutationState()){
                    newNode = pair.second;
                }else
                {
                    continue;
                }


            }else{
                if(newCn == cellCn){
                    newNode = pair.second;
                    break;
                }
            }
//



        }


        assert(newNode > -1);
        for(int i: phi.at(cellpair.first)){
            newPhi[newNode].insert(i);
        }


    }
    return newPhi;
}



/// given a genotype tree, update the genotypes for all nodes in the clonal tree
void ClonalTree::updateGenotype(int mut, Node u, GenotypeTree T)
{


    int mutState = T.mutationState();
    BoolNodeMap presence(_T, false);
    for(int state: D(_nodeToState[u]))
    {

        Node v =_stateToNode[state];
        CnaGenotype cn =_cnStates[v][_mut2seg[mut]];
        for(int i: T.D(mutState)){
            if (cn._x == T.x(i) && cn._y == T.y(i) && T.z(i) > 0)
            {
                presence[v] = true;
                _genotypes[v][mut] = Genotype(T.x(i), T.y(i), T.xbar(i), T.ybar(i));
                break;

            }
        }


    }
    for(NodeIt v(_T); v != lemon::INVALID; ++v)
    {
        if(!presence[v])
        {
            CnaGenotype cn =_cnStates[v][_mut2seg[mut]];
            _genotypes[v][mut] = Genotype(cn._x, cn._y, 0, 0);
        }
    }

}



/// given a genotype tree, update the genotypes for all nodes in the clonal tree
void ClonalTree::updateGenotype(IntSet muts, Node u, GenotypeTree T)
{


    if (!muts.empty())
    {

        //assume that all SNVs are in the same segment since they have the same genotype tree
        std::set<int>::iterator it = muts.begin();
        int element = *it;

        int mutState = T.mutationState();
        BoolNodeMap presence(_T, false);
        for (int state: D(_nodeToState[u]))
        {
            CnaGenotype cn = _cnStates[u][_mut2seg[element]];
            Node v = _stateToNode[state];
            presence[v] = false;
            for (int i: T.D(mutState))
            {
                if (cn._x == T.x(i) && cn._y == T.y(i) && T.z(i) > 0)
                {
                    presence[v] = true;
                    for (int mut: muts)
                    {
                        _genotypes[v][mut] = Genotype(T.x(i), T.y(i), T.xbar(i), T.ybar(i));

                    }
                    break;


                }
            }


        }

        for (NodeIt v(_T); v != lemon::INVALID; ++v)
        {
            if (!presence[v])
            {
                for(int mut: muts){
                    CnaGenotype cn = _cnStates[v][_mut2seg[mut]];
                    _genotypes[v][mut] = Genotype(cn._x, cn._y, 0, 0);
                }

            }
        }
    }


}

double ClonalTree::computeGenotypeLikelihood(int j,IntSetMap phi, const Data& D, double alpha, double lambda)
{
    double cnaCost = 0.0;
    double snvCost = 0.0;
    for(NodeIt u(_T); u != lemon::INVALID; ++u)
    {
        int nodeID = state(u);


        auto it = phi.find(nodeID);
        if (it != phi.end() && !it->second.empty()) {
            Genotype g = _genotypes[u][j];
            double vaf = g.vaf();
            double adjVaf = vaf * (1 - alpha) + (1 - vaf) * (alpha / 3);
            const IntSet& cells = it->second;
            for(int i: cells){
                cnaCost += abs(g._x - D._X[i][_mut2seg[j]])  + abs(g._y - D._Y[i][_mut2seg[j]]);
                snvCost += -1 * binomialProbability(D._Total[i][j],D._Var[i][j], adjVaf);
            }

        }

    }

    return snvCost + lambda * cnaCost;

}
GenotypeTree ClonalTree::getGenotypeTree(int snv)
{
    GenotypeTree::GenotypeEdgeSet edges;
    for(ArcIt e(_T); e != lemon::INVALID; ++e){
        Node parent = _T.source(e);
        Node child = _T.target(e);
        if(_genotypes[parent][snv] != _genotypes[child][snv])
        {
            edges.insert({_genotypes[parent][snv], _genotypes[child][snv]});
        }
    }

    return GenotypeTree(edges);

}

CnaTree ClonalTree::getCNATree(int segment)
{
    CnaTree::CnaEdgeSet edges;
    for(ArcIt e(_T); e != lemon::INVALID; ++e){
        Node parent = _T.source(e);
        Node child = _T.target(e);
        if(_cnStates[parent][segment] != _cnStates[child][segment])
        {
            edges.insert({_cnStates[parent][segment], _cnStates[child][segment]});
        }
    }

    return CnaTree(edges);

}

double ClonalTree::nodeCNAcost(const Node& u, const IntSet& cells,
                               const IntMatrix& xobs, const IntMatrix& yobs)
{
    double cost = 0.0;
    for(int ell: _segments)
    {
        CnaGenotype cn = _cnStates[u][ell];
        for(int i: cells){
            int x = xobs[i][ell];
            int y = yobs[i][ell];
            cost += abs(cn._x - x) + abs(cn._y - y);
        }
    }
    return cost;

}


double ClonalTree::nodeSNVcost(const Node& u, const IntSet& cells,
                               const IntMatrix& var, const IntMatrix& total,
                               double alpha)
{
    double cost = 0.0;
    for(int j: _muts)
    {
        Genotype g = _genotypes[u][j];
        double vaf = g.vaf();
        double adjVaf = vaf * (1 - alpha) + (1 - vaf) * (alpha / 3);
        assert(adjVaf >= 0);
        for(int i: cells)
        {
            int v = var[i][j];
            int t = total[i][j];
            if(t == 0){
                cost += 0;
            } else{

                cost += binomialProbability(t,v, adjVaf);

            }
        }
    }
    return cost;
}

std::pair<double, IntSetMap>  ClonalTree::assignCells(const IntSet& cells,
                                              const Data& D,
                                              double alpha,
                                              double lambda){
    IntMap phi;
    IntSetMap phiInv;
    double likelihood = 0.0;
    for (int i : cells) {
        double minCost = INFINITY;
        Node bestNode = lemon::INVALID;

        for (NodeIt u(_T); u != lemon::INVALID; ++u) {
            double snvCost = -1 * nodeSNVcost(u, {i}, D._Var, D._Total, alpha);
            double cnaCost = nodeCNAcost(u, {i}, D._X, D._Y);
            double nodeCost = snvCost + lambda * cnaCost;

            if (nodeCost < minCost) {
                bestNode = u;
                minCost = nodeCost;
            }
        }

        if (bestNode != lemon::INVALID) {
            likelihood += minCost;
            phi[i] = _nodeToState[bestNode];
        }
    }

    for (const auto& pair : phi) {
        phiInv[pair.second].insert(pair.first);
    }

    return {likelihood, phiInv};
}


double ClonalTree::binomialProbability(int n, int k, double p) {
    // Create a binomial distribution object
    boost::math::binomial_distribution<double> binom_dist(n, p);
    double prob = boost::math::pdf(binom_dist, k);
    return std::log(prob);
}

double ClonalTree::computeLikelihood(const IntSetMap& phi,
                                     const Data& D,
                                     double alpha,
                                     double lambda) {
    double cnaCost = 0.0;
    double snvCost = 0.0;

    for (NodeIt u(_T); u != lemon::INVALID; ++u)
    {
        int state = _nodeToState[u];
        auto it = phi.find(state);
        if (it != phi.end() && !it->second.empty()) {
            const IntSet& cells = it->second;
            cnaCost += nodeCNAcost(u, cells, D._X, D._Y);
            snvCost += -1 * nodeSNVcost(u, cells, D._Var, D._Total, alpha);
        }
    }

    return snvCost + lambda * cnaCost;
}

std::pair<double, IntSetMap> ClonalTree::optimize(const IntSet &cells, const Data& D,
                                                  double alpha, double lambda)
{

    std::pair<double, IntSetMap> cellAssignPair = assignCells(cells, D, alpha, lambda);
    double cost = computeLikelihood(cellAssignPair.second, D, alpha, lambda);
    double likelihood = cellAssignPair.first;
    double prevLikelihood;
//    IntSetMap  phi =

    do{
        prevLikelihood= likelihood;
        assignGenotypes(cellAssignPair.second, D, alpha, lambda);
        cost = computeLikelihood(cellAssignPair.second, D, alpha, lambda);
        for(NodeIt u(_T); u != lemon::INVALID; ++u){
            for(int j: _muts){
                Genotype g = _genotypes[u][0];
                std::cout << "Node " << state(u) << ":" << g._x << "," <<g._y << ","<< g._xbar <<"," << g._ybar << std::endl;
            }

        }
        cellAssignPair = assignCells(cells, D, alpha, lambda);

        likelihood = cellAssignPair.first;
        cost = computeLikelihood(cellAssignPair.second, D, alpha, lambda);



    }while(abs(prevLikelihood - likelihood) > 1);

    return cellAssignPair;
}


//std::ostream& operator<<(std::ostream& out, const ClonalTree& T)
//{
//    // output pi
////    out << -1;
////    for (int i = 1; i < T.k(); ++i)
////    {
////        out << " " << T.parent(i);
////    }
////    out << std::endl;
//
//    // output vertices
//    for (int i = 0; i < T.k(); ++i)
//    {
//        Node v_i = T.node(i);
//        out << "Node " << i << ":" << T._genotypes[v_i][0]._x
//            << " " <<T._genotypes[v_i][0]._y
//            << " " << T._genotypes[v_i][0]._xbar
//            << " " << T._genotypes[v_i][0]._ybar;
////        for (double s : T._s[v_i])
////        {
////            out << " " << s;
////        }
//        out << std::endl;
//    }
//
//    // output labels
//    out << T.label(0);
//    for (int i = 1; i < T.k(); ++i)
//    {
//        out << " " << T.label(i);
//    }
//    out << std::endl;
//
//    return out;
//}

int main() {
    std::vector<IntPair> edges = {{0, 1}, {0, 2}, {2, 3}};
    int k = 4;

    std::map<int, std::map<int, Genotype>> genotypes = {
            {0, {{0, Genotype {1, 1, 0, 0}}}},
            {1, {{0, Genotype {1, 1, 0, 0}}}},
            {2, {{0, Genotype {2, 1, 1, 0}}}},
            {3, {{0, Genotype {3, 1, 1, 0}}}},

    };

    int m = 2;
    for(int i= 0; i < k; ++i)
    {
        Genotype g = genotypes[i][0];
        for(int j = 1; j < m; ++j){
            genotypes[i][j] = Genotype(g._x, g._y, g._xbar, g._ybar);
        }

    }

    std::map<int, std::map<int, CnaGenotype>> cnStates = {
            {0, {{0, CnaGenotype {1, 1}}}},
            {1, {{0, CnaGenotype {1, 1}}}},
            {2, {{0, CnaGenotype {2, 1}}}},
            {3, {{0, CnaGenotype {3, 1}}}}
    };
    IntMap mut2seg = {{0, 0}, {1, 0}};
    IntSet mutClust = {2};
    ClonalTree tree(edges, genotypes, cnStates, mut2seg, mutClust);
    IntMatrix var(3, std::vector<int>(2));
    IntMatrix tot (3, std::vector<int>(2));
    IntMatrix xobs (3, std::vector<int>(1));
    IntMatrix yobs (3, std::vector<int>(1));
    var[0][0] = 0;
    var[1][0] = 20;
    var[2][0] = 20;

    tot[0][0] = 10;
    tot[1][0] = 40;
    tot[2][0] = 40;

    var[0][1] = 0;
    var[1][1] = 20;
    var[2][1] = 20;

    tot[0][1] = 10;
    tot[1][1] = 40;
    tot[2][1] = 40;

//    for (const auto& row : var) {
//        for (const auto& elem : row) {
//            std::cout << elem << " ";
//        }
//        std::cout << std::endl;
//    }

    xobs[0][0] = 1;
    xobs[1][0] = 2;
    xobs[2][0] = 3;

    yobs[0][0] = 1;
    yobs[1][0] = 1;
    yobs[2][0] = 1;
    IntMap  nodes = tree.getNodes();
    IntSetMap phi;

    phi[0] = {0};
    phi[2] = {2};
    phi[3] = {1};


    double alpha = 0.00;
    double lambda = 1000.0;
    IntSet cells = {0,1,2};
    Data D = Data(var, tot, xobs, yobs);
    double like = tree.computeLikelihood(phi, D,  alpha, lambda);
    std::cout << "Likelihood: " << like  << :: std::endl;

    std::pair<double, IntSetMap> cellAss= tree.optimize(cells, D, alpha, lambda);

//    std::pair<double, IntSetMap> cellAss = tree.assignCells(cells,D,0.001, 1000.0 );
    IntSetMap  pho = cellAss.second;
    std::cout << "Optimal likelihood: " << cellAss.first << std::endl;
    for(auto& pair: phi)
    {
        std::cout << "Node: " << pair.first <<  std::endl;
        for(int i: pair.second){
            std::cout << "cell: " << i << std::endl;
        }
    }
    return 0;
}