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
#include "singlecellgeneration.h"
#include <fstream>
#include <string>

int main(int argc, char** argv)
{
  int n = 1000;
  int m = 2;
  int k = 10;
  int kk = 5;
  int l = 5;
  int seed = 0;
  double expPurity = 0.8;
  double minProp = 0.05;
  int num_cells = 1000;
  double read_depth = .05;
  double alpha_fp = .001;
  std::string out_dir = "results";
  std::string dotFilename;
  std::string inputStateTreeFilename;
  std::string outputTreeFilename = "";
  std::string outputNodeFilename = "";
  std::string outputProportionFilename = "";
  bool removeUnsampledNodes = false;

  lemon2::ArgParser ap(argc, argv);
  ap.refOption("S", "Input CNA tree file", inputStateTreeFilename, true)
    .refOption("kk", "Number of truncal segments (default: 5)", kk, false)
    .refOption("s", "Random number generator seed (default: 0)", seed, false)
    .refOption("purity", "Expected purity (default: 0.8)", expPurity, false)
    .refOption("minProp", "Minimum desired clone proportion (default: 0.05)", minProp, false)
    .refOption("n", "Number of SNVs (default: 1000)", n, false)
    .refOption("m", "Number of samples (default: 2)", m, false)
    .refOption("k", "Number of segments (default: 10)", k, false)
    .refOption("l", "Number of mutation clusters (default: 5)", l, false)
    .refOption("STree", "Output filename for tree (default: none)", outputTreeFilename, false)
    .refOption("SNode", "Output filename for csv with information about nodes, segments, and mutaitons (default: none)", outputNodeFilename, false)
    .refOption("SProportions", "Output filename for csv with information about nodes, segments, and mutaitons (default: none)", outputProportionFilename, false)
    .refOption("dot", "Graphviz DOT output filename (default: '', no output)", dotFilename, false)
    .refOption("r", "Remove unsampled nodes", removeUnsampledNodes, false)
    .refOption("num_cells", "The number of cells to simulate with the single cell generation (default: 1000)", num_cells, false)
    .refOption("read_depth", "The read_depth for the single cell generation (default: .05)", read_depth, false)
    .refOption("alpha_fp", "The sequencing error for single cell generation (default .001)", alpha_fp, false)
    .refOption("out_dir", "The output directory for single cell generation (default: results)", out_dir, false);


  ap.parse();

  g_rng.seed(seed);

  if (!inputStateTreeFilename.empty())
  {
    std::ifstream inS(inputStateTreeFilename.c_str());
    if (!inS.good())
    {
      std::cerr << "Error: could not open '" << inputStateTreeFilename << "' for reading." << std::endl;
      return 1;
    }
    else
    {
      try
      {
        CnaGraph::readCnaTrees(inS);
      }
      catch (std::runtime_error& e)
      {
        std::cerr << e.what() << std::endl;
        return 1;
      }
    }
  }

  if (kk > k)
  {
    std::cerr << "Error: kk > k" << std::endl;
    return 1;
  }

  std::list<CnaTree> cnaTrees;
  for (int i = 0; i < kk;)
  {
    CnaTree T = CnaGraph::sampleCnaTree();
    if (T.truncal())
    {
      cnaTrees.push_back(T);
      ++i;
    }
  }

  for (int i = kk; i < k; ++i)
  {
    CnaTree T = CnaGraph::sampleCnaTree();
    cnaTrees.push_back(T);
  }

  Phylogeny phylo;
  for (const auto& T : cnaTrees)
  {
    phylo.addSegment(T, T.truncal());
  }

  phylo.sampleMutations(n, l); 
  phylo.sampleProportions(m, expPurity, minProp);


/*   phylo.writeDOT(std::cout); 
  phylo.writeTree(std::cout, outputTreeFilename);
  phylo.writeNodeFile(std::cout, outputNodeFilename);
  phylo.writeProportionFile(std::cout, outputProportionFilename, m); */


  if (removeUnsampledNodes) 
  {
    Phylogeny newPhylo = phylo.removeUnsampledNodes();
    //newPhylo.writeDOT(std::cout); //this line
    newPhylo.writeTree(std::cout, outputTreeFilename);
    newPhylo.writeNodeFile(std::cout, outputNodeFilename);
    newPhylo.writeProportionFile(std::cout, outputProportionFilename, m);
/*     std::cout << newPhylo; 
    if (!dotFilename.empty())
    {
      std::ofstream outDOT(dotFilename);
      newPhylo.writeDOT(outDOT);
      outDOT.close();
    } */
  }
/*   else
  {
    std::cout << phylo;
    if (!dotFilename.empty())
    {
      std::ofstream outDOT(dotFilename);
      phylo.writeDOT(outDOT);
      outDOT.close();
    }
  } */
  

  
  for (int i = 0; i < m; i++) {
    SingleCell sc(num_cells, read_depth, alpha_fp, out_dir, k, m); 
    sc.loadData(std::cout, outputProportionFilename, outputNodeFilename);
    sc.generateECDF(std::cout, i);
    sc.initializeSCS(std::cout); 
    sc.generateCells(std::cout, i, g_rng);
    sc.printSCS(std::cout, i);
  }

  return 0;
}