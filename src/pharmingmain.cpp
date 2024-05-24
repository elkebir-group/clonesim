//
// Created by Leah Weber on 5/20/24.
//

#include "utils.h"
#include "clonaltree.h"
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
    phi[2] = {1};
    phi[3] = {2};


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