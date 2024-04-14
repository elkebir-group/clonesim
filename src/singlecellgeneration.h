#ifndef SINGLECELLGENERATION_H
#define SINGLECELLGENERATION_H

#include "utils.h"
#include <map>

class SingleCell {
public:
    //Default constructor
    SingleCell(int numSCS, double read_depth, double alpha_fp, std::string outdir, int numSegments, int numSamples,
               double cna_error);

    void main(std::ostream &out, std::string &input_file_dir, int sample, double cnaErrorRate, double plsOne, double minusOne);

    void loadData(std::ostream &out, std::string &input_files, double cnaErrorRate, double plsOne, double minusOne);

    void addNoiseCNA(double cnaErrorRate, double plsOne, double minusOne);

    int gaussianDraw(int mean, double errorRate);

    void generateECDF(int i);

    void initializeSCS();

    void generateCells(int sample);

    int sampleSingleCells(int sample);

    std::pair<int, int> draw(int mut_alleles, int total_cn);

    int binomialdraw(double p, int n);

    void printSCS(std::ostream &out, int sample);


private:
    const int _NUMSCS;
    double READ_DEPTH;
    int _numSegments;
    double _cnaError;
    //_seed;
    std::vector<std::vector<int>> _varReads; //The number of variant reads at a mutation locus in each single cell. Dimensions: #cells x #mutations
    std::vector<std::vector<int>> _refReads; // The number of non-mutated reads at a mutation locus in each single cell. Dimensions: #cells x #mutations
    std::vector<std::vector<int>> _totReads; // The number of total reads at a mutation locus in each single cell. DImensions: #cells x #mutations
    std::vector<int> _cellLabels; //the clone ID that each cell belongs to (dimension: # of cells)
    std::vector<std::vector<double>> _SCS_PREV; //the proportion of each node in each sample (dimensions: # of nodes x # of samples)
    std::vector<std::vector<int>> _NODE_INFORMATION;
    std::vector<std::vector<double>> _ecdf; //cdf of clonal proportions (dimensions: 1 x #clones)
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
    double _READ_DEPTH;
    double _alpha_fp;


};


#endif //SINGLECELLGENERATION_H