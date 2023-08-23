//
// Created by Mohammed El-Kebir on 8/13/23.
//

#include "utils.h"
#include "cnagraph.h"
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
  int maxTotalCN = 4;
  int minTotalCN = 0;
  int rootX = 1;
  int rootY = 1;
  std::string inputCnaTreeFilename;
  std::string outputCnaTreeFilename;

  lemon2::ArgParser ap(argc, argv);
  ap.refOption("maxXY", "Maximum number of maternal/paternal copies (default: 2)", maxXY, false)
    .refOption("minTotalCN", "Minimum total copy number (default: 0)", minTotalCN, false)
    .refOption("maxTotalCN", "Maximum total copy number (default: 4)", maxTotalCN, false)
    .refOption("maxCN", "Maximum number of copy number events (default: 2)", maxCN, false)
    .refOption("S", "Input state tree file", inputCnaTreeFilename, false)
    .refOption("SS", "Output state tree file", outputCnaTreeFilename, false)
    .refOption("rootX", "Number of maternal copies in the root (default: 1)", rootX, false)
    .refOption("rootY", "Number of paternal copies in the root (default: 1)", rootY, false);
  ap.parse();

  if (!inputCnaTreeFilename.empty())
  {
    std::ifstream inS(inputCnaTreeFilename.c_str());
    if (!inS.good())
    {
      std::cerr << "Error: could not open '" << inputCnaTreeFilename << "' for reading." << std::endl;
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

          maxCopy = std::max(maxCopy, std::max(xy.first, xy.second));
        }
        std::cerr << "..." << std::flush;

        CnaGraph::getCnaTrees(L, rootX, rootY);

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
        if (minTotalCN <= x + y &&  x + y <= maxTotalCN)
        {
          cnStates.push_back(IntPair(x, y));
          if (x != y)
          {
            cnStates.push_back(IntPair(y, x));
          }
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
        maxCopy = std::max(maxCopy,
                           std::max(cnStates[enumState[i]].first,
                                    cnStates[enumState[i]].second));
      }

      for (int i = 0; i < maxCN; ++i)
      {
        if (i > 0)
          std::cerr << " ";
        const IntPair& xy = cnStates[enumState[i]];
        std::cerr << "(" << xy.first << "," << xy.second << ")";
      }
      std::cerr << "..." << std::flush;

      std::cerr << " Enumerated " << CnaGraph::getCnaTrees(L, rootX, rootY).size()
                << " CNA trees." << std::endl;
    } while (next(maxXY, maxCN, cnStates.size(), minEnumState, maxEnumState, enumState));
  }

  if (outputCnaTreeFilename.empty())
  {
    CnaGraph::writeCnaTrees(std::cout);
  }
  else
  {
    std::ofstream out(outputCnaTreeFilename.c_str());
    CnaGraph::writeCnaTrees(out);
    out.close();
  }

  return 0;
}
