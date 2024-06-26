/*
 * simulatemain.cpp
 *
 *  Created on: 13-aug-2023
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "cnagraph.h"
#include "arg_parser.h"
#include "phylogeny.h"
#include <fstream>
#include <string>

int main(int argc, char **argv) {
    int n = 1000;
    int m = 2;
    int k = 10;
    int kk = 5;
    int l = 5;
    int seed = 0;
    int num_tries = 1000;
    double dirich_param = 2;
    double expPurity = 0.8;
    double minProp = 0.05;
    bool _f = false;
    bool restrictLoss = false;
    std::string out_dir = "results";
    std::string dotFilename;
    std::string inputStateTreeFilename;
    std::string _output_file_dir = "";
    bool removeUnsampledNodes = false;
    bool uniform = false;
    double threshold = 0.05;

    lemon2::ArgParser ap(argc, argv);
    ap.refOption("S", "Input CNA tree file", inputStateTreeFilename, true)
            .refOption("dirich_param", "symmetric concentration parameter for Dirichlet distribution (default 2)", dirich_param, false)
            .refOption("kk", "Number of truncal segments (default: 5)", kk, false)
            .refOption("s", "Random number generator seed (default: 0)", seed, false)
            .refOption("purity", "Expected purity (default: 0.8)", expPurity, false)
            .refOption("minProp", "Minimum desired clone proportion (default: 0.05)", minProp, false)
            .refOption("threshold", "Minimum threshold for SNV proportions (default: 0.05)", threshold, false)
            .refOption("n", "Number of SNVs (default: 1000)", n, false)
            .refOption("m", "Number of samples (default: 2)", m, false)
            .refOption("k", "Number of segments (default: 10)", k, false)
            .refOption("l", "Number of mutation clusters (default: 5)", l, false)
            .refOption("dot", "Graphviz DOT output filename (default: '', no output)", dotFilename, false)
            .refOption("r", "Remove unsampled nodes", removeUnsampledNodes, false)
            .refOption("f", "Whether to output files", _f, false)
            .refOption("output_file_dir", "The directory for where to write output files", _output_file_dir, false)
            .refOption("num_tries", "The number of tries for sampling mutation rejection sampling (default 1000)", num_tries, false)
            .refOption("restrictLoss", "Whether to restrict copy number loss (default false)", restrictLoss, false)
            .refOption("uniform", "use uniform distribution for mutation assignments", uniform, false);
            //.refOption("dirichletParam", "The parameter for the dirichlet (default 1)", dirich_param, false)



    ap.parse();

    g_rng.seed(seed);

    if (!inputStateTreeFilename.empty()) {
        std::ifstream inS(inputStateTreeFilename.c_str());
        if (!inS.good()) {
            std::cerr << "Error: could not open '" << inputStateTreeFilename << "' for reading." << std::endl;
            return 1;
        } else {
            try {
                CnaGraph::readCnaTrees(inS);
            }
            catch (std::runtime_error &e) {
                std::cerr << e.what() << std::endl;
                return 1;
            }
        }
    }

    std::cerr << "CNA trees read in, constructing clonal tree..." << std::endl;

    if (kk > k) {
        std::cerr << "Error: kk > k" << std::endl;
        return 1;
    }

    try {
        std::list<CnaTree> cnaTrees;
        for (int i = 0; i < kk;) {
            CnaTree T = CnaGraph::sampleCnaTree();
            if (T.truncal() && (!restrictLoss | !T.hasLoss())) {
                cnaTrees.push_back(T);
                ++i;
            }
        }

        for (int i = kk; i < k;) {
            CnaTree T = CnaGraph::sampleCnaTree();
            if (!restrictLoss | !T.hasLoss()) {
                cnaTrees.push_back(T);
                ++i;
            }
        }

        Phylogeny phylo;
        for (const auto &T: cnaTrees) {
            phylo.addSegment(T, T.truncal());
        }
        phylo.createIndex();
        //std::ofstream outTree("/Users/annahart/CLionProjects/clonesim/build/test/tree.txt");
        //phylo.writeTree(outTree);

        phylo.sampleMutations(n, l, num_tries, dirich_param, uniform, threshold);


        std::cerr << "Clonal tree constructed, sampling proportions..." << std::endl;

        phylo.sampleProportions(m, expPurity, minProp);

        std::cerr << "Finished sampling proportions";

        std::cerr << "Removing unsampled nodes..." << std::endl;
        if (removeUnsampledNodes) {
            Phylogeny newPhylo = phylo.removeUnsampledNodes();
            std::cout << newPhylo;

            if (_f) {
                std::cerr << "Writing clonal tree output files..." << std::endl;
                newPhylo.writeFiles(std::cout, _output_file_dir, m);
            }

            if (!dotFilename.empty()) {
                std::ofstream outDOT(dotFilename);
                newPhylo.writeDOT(outDOT);
                outDOT.close();
            }
        } else {

            phylo.writeDOT(std::cout);

            if (_f) {
                std::cerr << "Writing clonal tree output files..." << std::endl;
                phylo.writeFiles(std::cout, _output_file_dir, m);
            }


            if (!dotFilename.empty()) {
                std::ofstream outDOT(dotFilename);
                phylo.writeDOT(outDOT);
                outDOT.close();
            }
        }


    }
    catch (std::runtime_error &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}