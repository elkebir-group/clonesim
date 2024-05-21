#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include "gencnatrees.h"
#include "../cnagraph.h"
#include "../cnatree.h"
#include "../clonaltree.h"

namespace py = pybind11;

typedef std::tuple<int, int, int, int> genotype;
typedef std::map<int, genotype> snvGenotype;
typedef std::map<int, snvGenotype> genotypes;


//wrapper to convert CnaTree:CnaEdgeSetSet to python friendly format (python doesn't seem to like sets)
const std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> genCNATrees(const IntPairSet& L,
                                                    int root_x, int root_y){

    CnaTree::CnaEdgeSetSet scriptS = CnaGraph::getCnaTrees(L, root_x, root_y);
    std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>> scriptT;

    for (auto tree : scriptS) {
        std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> tr;

        for (auto edge : tree) {
            // Extract parent and child values
            int parent_x = edge.first._x;
            int parent_y = edge.first._y;
            int child_x = edge.second._x;
            int child_y = edge.second._y;

     
            // // Create pairs of integers for parent and child
            std::pair<int, int> par(parent_x, parent_y);
            std::pair<int, int> chi(child_x, child_y);

            // Create a pair of pairs for the edge
            std::pair<std::pair<int, int>, std::pair<int, int>> edge_pair(par, chi);

            // Insert the pair of pairs into the vector
            tr.push_back(edge_pair);
        }
  
        // Insert the vector of pairs of pairs into the outer vector
        scriptT.push_back(tr);
    }


    return scriptT;

 }
//  

 std::vector<std::vector<std::pair<std::tuple<int, int,int, int>, std::tuple<int, int, int, int>>>> genGenoTrees(std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> cnatree){
        
        CnaTree::CnaEdgeSet tree;
//        std::cout << "Number of edges: " << cnatree.size() << std::endl;

        for(auto edge: cnatree){
            std::pair<int, int> par = edge.first; 
            std::pair<int, int> child = edge.second; 
//            std::cout << par.first << "," << par.second << "->" << child.first << "," << child.second <<std::endl;
            CnaTree::CnaGenotype pGeno(par.first, par.second);
            CnaTree::CnaGenotype cGeno(child.first, child.second);
            CnaTree:: CnaEdge ed = std::make_pair(pGeno, cGeno);
            tree.insert(ed);
        }

        CnaTree cntree(tree);


        

//	std::cout << "About to enumerate..." << std::endl;

        GenotypeTree::GenotypeEdgeSetSet result;
        cntree.enumerateGenotypeTrees(result);
 


         std::vector<std::vector<std::pair<std::tuple<int, int,int, int>, std::tuple<int, int, int, int>>>> scriptT;

        for(auto gtree: result){
            std::vector<std::pair<std::tuple<int, int,int, int>, std::tuple<int, int, int, int>>> geno_tree;
            for(auto e: gtree){
                GenotypeTree::Genotype parent= e.first;
                GenotypeTree::Genotype child= e.second;
                geno_tree.push_back(std::make_pair( 
                    std::make_tuple(parent._x, parent._y, parent._xbar, parent._ybar),
                    std::make_tuple(child._x, child._y, child._xbar, child._ybar) 
                    ));


            }
            scriptT.push_back(geno_tree);

        }
   
//	std::cout << "Done" << std::endl;
        

        return scriptT;


 }

 std::tuple<double, std::map<int,int> , genotypes> optimizeClonalTree( std::vector<std::pair<int,int>> edges,
                                                               genotypes genos,
                                                               std::map<int, std::map<int, std::pair<int,int>>> cnstates,
                                                        IntMap mut2seg,
                                                        IntSet mutClusts,
                                                        IntSet cells,
                                                        IntMatrix var,
                                                        IntMatrix total,
                                                        IntMatrix xobs,
                                                        IntMatrix yobs,
                                                        double alpha,
                                                        double lambda){

     std::map<int,int> phi;
     phi[0] =1;

    std::map<int, IntToGenotypeMap> inGenos;
//
    for(auto nodeGeno: genos)
    {
        int node = nodeGeno.first;
        for(auto mutGeno: nodeGeno.second)
        {
            std::tuple<int,int,int,int> g = mutGeno.second;
            inGenos[node][mutGeno.first] =Genotype(std::get<0>(g),std::get<1>(g),
                                                   std::get<2>(g),std::get<3>(g));
        }

    }

//     for(auto pair: inGenos){
//         int n = pair.first;
//         for(auto mutG: pair.second)
//         {
//             Genotype g = inGenos[n][mutG.first];
//             std::cout << "Node " << n << ":" << g._x << "," <<g._y << ","<< g._xbar <<"," << g._ybar << std::endl;
//         }
//
//     }

//
    std::map<int, IntToCnaGenotypeMap> inCnStates;
    for(auto nodeCn: cnstates)
    {
        int node = nodeCn.first;
        for(auto segCn: nodeCn.second)
        {
            inCnStates[node][segCn.first] = CnaGenotype(segCn.second.first,segCn.second.second );
        }

    }
     Data D(var, total, xobs, yobs);

    ClonalTree ct(edges, inGenos, inCnStates, mut2seg, mutClusts);


    std::pair<double,IntSetMap> result =  ct.optimize(cells, D, alpha, lambda);


     for(auto pair: result.second)
     {
         int n = pair.first;
         for(int i: pair.second){
             phi[i] =n;
         }
     }

    genotypes G = ct.getGenotypes();

    IntSetMap phiInv = result.second;
//    return std::make_tuple(100.00, phi, genos);



    return std::make_tuple(result.first, phi, G);

}



PYBIND11_MODULE(clonelib, m) {

    m.def("get_cna_trees", &genCNATrees, "Get CNA trees as a Python list of lists of tuples of int pair tuples (tree edges)");
    m.def("get_genotype_trees", &genGenoTrees, "Get genotype trees as a Python list of lists of tuples of int 4 tuples (tree edges)");
    m.def("optimize_clonal_tree", &optimizeClonalTree, "Optimize a clonal tree.");

}
