//
// Created by Mohammed El-Kebir on 9/3/23.
//

#include "utils.h"
#include "cnagraph.h"
#include "arg_parser.h"
#include "phylogeny.h"
#include <fstream>

int main(int argc, char** argv)
{
  if (argc != 2 && argc != 1)
  {
    std::cerr << "Usage: " << argv[0] << " (<phylo>)" << std::endl;
    return 1;
  }

  Phylogeny T;

  try
  {
    if (argc == 1)
    {
      std::cin >> T;
    } else
    {
      assert(argc == 2);
      std::ifstream inPhylo(argv[1]);
      if (!inPhylo.good())
      {
        std::cerr << "Could not open '" << argv[1] << "' for reading" << std::endl;
        return 1;
      }
      inPhylo >> T;
      inPhylo.close();
    }
  }
  catch (std::runtime_error& e)
  {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  T.writeDOT(std::cout);

  return 0;
}