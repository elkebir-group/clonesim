/*
 * mergegenotypetreesmain.cpp
 *
 *  Created on: 10-jun-2019
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "genotypegraph.h"
#include <lemon/arg_parser.h>
#include <fstream>

int main(int argc, char** argv)
{
  if (argc < 2)
  {
    std::cerr << "Usage: " << argv[0] << " <genotype_tree_file_1> ... <genotype_tree_file_n>" << std::endl;
    return 0;
  }

  for (int i = 1; i < argc; ++i)
  {
    std::ifstream inS(argv[i]);
    if (inS.good())
    {
      try
      {
        GenotypeGraph::readGenotypeTrees(inS);
      }
      catch (std::runtime_error& e)
      {
        std::cerr << e.what() << std::endl;
        return 1;
      }
    }
    else
    {
      std::cerr << "Error: Could not open '" << argv[i] << "' for reading." << std::endl;
      return 1;
    }
  }

  GenotypeGraph::writeGenotypeTrees(std::cout);
  
  return 0;
}
