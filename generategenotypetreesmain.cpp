/*
 * generategenotypetreesmain.cpp
 *
 *  Created on: 28-mar-2019
 *      Author: M. El-Kebir
 */

#include "utils.h"
#include "genotypegraph.h"
#include "cnagraph.h"
#include "cnatree.h"
#include "arg_parser.h"
#include <fstream>

std::vector<IntPairSet> parseCopyStateInput(std::istream& in)
{
  std::vector<IntPairSet> res;
  
  while (in.good())
  {
    std::string line;
    getline(in, line);
    
    if (line.empty()) continue;
    
    StringVector s;
    boost::split(s, line, boost::is_any_of(";"));
    
    IntPairSet L;
    for (const std::string& str : s)
    {
      StringVector ss;
      boost::split(ss, str, boost::is_any_of(","));
      
      if (ss.size() != 2)
      {
        throw std::runtime_error(getLineNumber() + "Error: Invalid copy-number state '" + str + "'");
      }
      
      try
      {
        int x = boost::lexical_cast<int>(ss[0]);
        int y = boost::lexical_cast<int>(ss[1]);
      
      L.insert(IntPair(x, y));
      }
      catch (boost::exception& e)
      {
        throw std::runtime_error(getLineNumber() + "Error: Invalid copy-number state '" + str + "'");
      }
    }
    res.push_back(L);
  }
  
  return res;
}

bool next(const int maxXY,
          const int maxCN,
          const int nrCnStates,
          const IntVector& minEnumState,
          const IntVector& maxEnumState,
          IntVector& enumState)
{
  assert(enumState.size() == maxCN);
  
  IntVector newEnumState = enumState;
  
  int i = maxCN - 1;
  for (; i >= 0; --i)
  {
    if (newEnumState[i] < maxEnumState[i] - 1)
    {
      ++newEnumState[i];
      break;
    }
  }
  
  if (i < 0)
  {
    return false;
  }
  else
  {
    for (int j = i + 1; j < maxCN; ++j)
    {
      newEnumState[j] = newEnumState[j-1] + 1;
    }
  }
  
  enumState = newEnumState;
  
  return true;
}

int main(int argc, char** argv)
{
  int maxXY = 2;
  int maxCN = 2;
  int rootX = 1;
  int rootY = 1;
  std::string inputStateTreeFilename;
  std::string outputStateTreeFilename;
  bool viaCnaTrees = true;
  
  lemon2::ArgParser ap(argc, argv);
  ap.refOption("maxXY", "Maximum number of maternal/paternal copies (default: 2)", maxXY, false)
    .refOption("maxCN", "Maximum number of copy number events (default: 2)", maxCN, false)
    .refOption("S", "Input genotype tree file", inputStateTreeFilename, false)
    .refOption("SS", "Output genotype tree file", outputStateTreeFilename, false)
//    .refOption("cna", "Use CNA trees as an intermediate", viaCnaTrees, false)
    .refOption("rootX", "Number of maternal copies in the root (default: 1)", rootX, false)
    .refOption("rootY", "Number of paternal copies in the root (default: 1)", rootY, false);
  ap.parse();
  
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
        GenotypeGraph::readGenotypeTrees(inS);
      }
      catch (std::runtime_error& e)
      {
        std::cerr << e.what() << std::endl;
        return 1;
      }
    }
  }
  
  if (!ap.files().empty())
  {
    for (const std::string& filename : ap.files())
    {
      std::ifstream in(filename.c_str());
      if (!in.good())
      {
        std::cerr << "Error: could not open '" << filename << "' for reading" << std::endl;
        return 1;
      }
      
      std::vector<IntPairSet> cnStates;
      try
      {
        cnStates = parseCopyStateInput(in);
      }
      catch (std::runtime_error& e)
      {
        std::cerr << e.what() << std::endl;
        return 1;
      }
      in.close();
      
      for (const IntPairSet& L : cnStates)
      {
        int maxCopy = 0;
        bool first = true;
        for (const IntPair& xy : L)
        {
          if (!first)
            std::cerr << " ";
          else
            first = false;
          
          std::cerr << "(" << xy.first << "," << xy.second << ")";
        }
        std::cerr << "..." << std::flush;
        
        GenotypeGraph::getStateTrees(L, rootX, rootY);
        
        std::cerr << " Done." << std::endl;
      }
      
      std::cerr << "Processed " << filename << "..." << std::endl;
    }
  }
  else
  {
    std::vector<IntPair> cnStates;
    for (int x = 0; x <= maxXY; ++x)
    {
      for (int y = x; y <= maxXY; ++y)
      {
        cnStates.push_back(IntPair(x, y));
        if (x != y)
        {
          cnStates.push_back(IntPair(y, x));
        }
      }
    }
    
    IntVector enumState(maxCN);
    for (int i = 0; i < maxCN; ++i)
    {
      enumState[i] = i;
    }
    
    IntVector minEnumState = enumState;
    
    IntVector maxEnumState(maxCN);
    for (int i = 0; i < maxCN; ++i)
    {
      maxEnumState[i] = cnStates.size() - (maxCN - i - 1);
    }
    
    do
    {
      IntPairSet L;
      int maxCopy = 0;
      for (int i = 0; i < maxCN; ++i)
      {
        L.insert(cnStates[enumState[i]]);
      }
      
      for (int i = 0; i < maxCN; ++i)
      {
        if (i > 0)
          std::cerr << " ";
        const IntPair& xy = cnStates[enumState[i]];
        std::cerr << "(" << xy.first << "," << xy.second << ")";
      }
      std::cerr << "..." << std::flush;

      if (viaCnaTrees)
      {
        const auto& cnaTrees = CnaGraph::getCnaTrees(L, rootX, rootY);
        for (const auto& cnaTree : cnaTrees)
        {
          CnaTree T(cnaTree);
          T.enumerateGenotypeTrees(GenotypeGraph::_dict[L]);
        }
      }
      else
      {
        if (viaCnaTrees)
        {
          const auto& cnaTrees = CnaGraph::getCnaTrees(L, rootX, rootY);
          for (const auto& cnaTree : cnaTrees)
          {
            CnaTree T(cnaTree);
            T.enumerateGenotypeTrees(GenotypeGraph::_dict[L]);
          }
        }
        else
        {
          GenotypeGraph::getStateTrees(L, rootX, rootY);
        }
      }
      
      std::cerr << " Done." << std::endl;
    } while (next(maxXY, maxCN, cnStates.size(), minEnumState, maxEnumState, enumState));
  }
  
  if (outputStateTreeFilename.empty())
  {
    GenotypeGraph::writeGenotypeTrees(std::cout);
  }
  else
  {
    std::ofstream out(outputStateTreeFilename.c_str());
    GenotypeGraph::writeGenotypeTrees(out);
    out.close();
  }
  
  return 0;
}
