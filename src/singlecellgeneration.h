#ifndef SINGLECELLGENERATION_H
#define SINGLECELLGENERATION_H

#include "utils.h"
#include <map>

class SingleCell {
public: 
    //Default constructor
    SingleCell(int numSCS, float read_depth, float alpha_fp, std::string outdir, int numSegments, int numSamples, std::mt19937 &gen, double cna_error);
    void loadData(std::ostream& out, std::string& proportionsFileName, std::string& outputNodeFilename) ;
    void addNoiseCNA();
    int gaussianDraw(int mean, double errorRate);
    void generateECDF(std::ostream& out, int i);
    void resizeArrays(std::ostream& out);
    void initializeSCS(std::ostream& out);
    void generateCells(std::ostream& out, int sample);
    int sampleSingleCells(std::ostream& out, int sample);
    std::pair<int, int> draw(std::ostream& out, int mut_alleles, int ref_alleles);
    int binomialdraw(float p, int n);
    void printSCS(std::ostream& out, int sample);
    void write_csv(std::string filename, std::vector<int> cell_ids,
                         std::vector<std::string> colnames,
                         std::vector<std::vector<int> > dataset, std::string delim);


private: 
     const int _NUMSCS;
     float READ_DEPTH;
     int _numSegments;
     double _cnaError;
     std::mt19937 _gen;
     //_seed;
     std::vector<std::vector<int>> _varReads; //The number of variant reads at a mutation locus in each single cell. Dimensions: #cells x #mutations 
     std::vector<std::vector<int>> _refReads; // The number of non-mutated reads at a mutation locus in each single cell. Dimensions: #cells x #mutations
     std::vector<std::vector<int>> _totReads; // The number of total reads at a mutation locus in each single cell. DImensions: #cells x #mutations
     std::vector<int> _cellLabels; //the clone ID that each cell belongs to (dimension: # of cells)
     std::vector<std::vector<float>> _SCS_PREV; //the proportion of each node in each sample (dimensions: # of nodes x # of samples)
     std::vector<std::vector<int>> _NODE_INFORMATION; 
     std::vector<std::vector<float>> _ecdf; //cdf of clonal proportions (dimensions: 1 x #clones)
     std::vector<std::vector<std::vector<int>>> _mutAlleles; // Gives mutated copies of allele @ the locus of mutation in that clone. Dimensions: #clones x # muatations x 2
     std::vector<std::vector<std::vector<int>>> _refAlleles; // Gives non-mutated copies of allele @ the locus of mutation in that clone. Dimensions: #clones x #mutations x 2
     std::vector<std::vector<std::vector<int>>> _totalAlleles; //Gives total copies of allele @ the locaus of mutation in that clone. Dimensions: #clones x #mutations x 2
     std::vector<std::vector<std::vector<int>>> _segmentCopyNumbers; //The copy number of each segemnt in each clone for x and y. Dimensions: #clones x #segments x 2
     std::vector<std::vector<int>> _segments; //the segment of each mutation in each clone. Dimensions: #clones x #mutations
     std::string outdir;
     std::map<std::string, std::string> _outFiles;
     int _numClones;
     int _numSamples;
     int _nodeInformationColumns;
     int _nodeInformationRows;
     int _numMutations;
     int _nodeCol;
     int _segmentCol;
     int _xCol;
     int _yCol;
     int _mCol;
     int _xbarCol;
     int _ybarCol;
     int _clusterIDCol;
     float _READ_DEPTH;
     float _alpha_fp;
     
  


    
};


#endif //SINGLECELLGENERATION_H