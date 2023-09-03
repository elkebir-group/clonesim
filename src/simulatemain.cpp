/*
 * simulatemain.cpp
 *
 *  Created on: 13-aug-2023
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "genotypetree.h"
#include "genotypegraph.h"
#include "cnagraph.h"
#include "arg_parser.h"
#include "phylogeny.h"
#include <fstream>

int main(int argc, char** argv)
{
  int n = 1000;
  int m = 2;
  int k = 10;
  int kk = 5;
  int l = 5;
  int seed = 0;
  std::string inputStateTreeFilename;
  std::string outputTreeFilename = "";
  std::string outputNodeFilename = "";

  lemon2::ArgParser ap(argc, argv);
  ap.refOption("S", "Input CNA tree file", inputStateTreeFilename, true)
    .refOption("kk", "Number of truncal segments (default: 5)", kk, false)
    .refOption("s", "Random number generator seed (default: 0)", seed, false)
    .refOption("n", "Number of SNVs (default: 1000)", n, false)
    .refOption("m", "Number of samples (default: 2)", m, false)
    .refOption("k", "Number of segments (default: 10)", k, false)
    .refOption("l", "Number of mutation clusters (default: 5)", l, false)
    .refOption("STree", "Output filename for tree (default: none)", outputTreeFilename, false)
    .refOption("SNode", "Output filename for csv with information about nodes, segments, and mutaitons (default: none)", outputNodeFilename, false);
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

  phylo.writeDOT(std::cout);
  phylo.writeTree(std::cout, outputTreeFilename);
  phylo.writeNodeFile(std::cout, outputNodeFilename);

  return 0;
}