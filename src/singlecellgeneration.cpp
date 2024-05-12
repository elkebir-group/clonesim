#include "singlecellgeneration.h"
#include <fstream>
#include <string>
#include <lemon/connectivity.h>
#include "utils.h"
#include <cmath>
#include <unordered_map>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/binomial_distribution.hpp>

SingleCell::SingleCell(int numSCS, double read_depth, double alpha_fp, std::string outdir, int numSegments,
                       int numSamples, double cna_error)
        : _NUMSCS(numSCS) //** HYPERPARAMETER
        , _READ_DEPTH(read_depth) //** HYPERPARAMETER
        , _numSegments(numSegments), _numSamples(numSamples), _cnaError(cna_error), outdir(outdir), _outFiles(
                {{"SPARSE",          outdir + "/sparse"},
                 {"CELLS",           outdir + "/cells"},
                 {"CELLASSIGNMENTS", outdir + "/cellAssignments"}}), _alpha_fp(alpha_fp)
        //, _seed
        , _varReads(), _refReads(), _totReads(), _cnCalls(numSCS, std::vector<std::pair<int, int>>(numSegments)), _cellLabels(), _SCS_PREV() //Note: the first column is the node ID
        , _NODE_INFORMATION(), _ecdf() //cdf of clonal proportions
        , _segmentCopyNumbers(), _segments(), _numClones(0), _nodeInformationRows(0), _nodeInformationColumns(0),
          _numMutations(0), _nodeCol(0), _segmentCol(1), _xCol(2), _yCol(3), _mCol(4), _xbarCol(5), _ybarCol(6),
          _clusterIDCol(7) {

}

void SingleCell::main(std::ostream &out, std::string &input_file_dir, int i, double cnaError, double plsOne) {
    std::cerr << "loading data" << std::endl;
    loadData(std::cout, input_file_dir);
    std::cerr << "generating ecdf" << std::endl;
    generateECDF(i);
    std::cerr << "initializing ecdf" << std::endl;
    initializeSCS();
    std::cout << "Done intializing." << std::flush;
    std::cerr << "generate cells" << std::endl;
    generateCells(i, cnaError, plsOne);
    std::cerr << "saving data" << std::endl;
    printSCS(std::cout, i);
    std::cerr << "sample  " << i << " complete!" << std::endl;
}

void SingleCell::loadData(std::ostream &out, std::string &input_file_dir) {
    int i = 0;
    int j = 0;
    std::string proportion_fn = input_file_dir + "/proportions.tsv";
    std::cout << proportion_fn << std::endl;
    std::string node_fn = input_file_dir + "/node.tsv";
    std::cout << node_fn << std::endl;
    std::ifstream file(proportion_fn);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open node probability file");
    }
    std::string line;
    getline(file, line); //discard header 
    std::string cell;
    double proportion = 0.0;
    while (std::getline(file, line)) {
        j = 0;
        std::vector<double> row;
        std::stringstream linestream(line);
        while (std::getline(linestream, cell, '\t')) {
            proportion = std::stof(cell);
            row.push_back(proportion);
            j++;
        }
        _SCS_PREV.push_back(row);
        i++;
    }
    _numClones = i;
    _numSamples = j - 1;
    file.close();

    ///////////////////////////

    int ii = 0;
    int jj = 0;

    std::ifstream file2(node_fn);
    if (!file2.is_open()) {
        throw std::runtime_error("Could not open node file");
    }
    std::string line2;
    getline(file2, line2); //discard header 
    std::string cell2;
    int value = 0;
    while (std::getline(file2, line2)) {
        jj = 0;
        std::vector<int> row2;
        std::stringstream linestream(line2);
        while (std::getline(linestream, cell2, '\t')) {
            value = std::stoi(cell2);
            row2.push_back(value);
            jj++;
        }
        _NODE_INFORMATION.push_back(row2);
        ii++;
    }
    _nodeInformationColumns = jj;
    _nodeInformationRows = ii;


/*     for (int aa = 0; aa < _nodeInformationRows; aa++) {
        for (int bb = 0; bb < _nodeInformationColumns; bb++) {
            out << _NODE_INFORMATION[aa][bb] << "\t";
        }
        out << "\n";
    } 

    out << "\n\n*******\n\n"; */


}


void SingleCell::addCNANoise(int cell, int segment, double cnaErrorRate, double plsOne) {
    int clone = _cellLabels[cell];
    int xcopy = _segmentCopyNumbers[clone][segment][0];
    int ycopy = _segmentCopyNumbers[clone][segment][1];
    boost::random::uniform_01<> myrand;
    if (myrand(g_rng) < cnaErrorRate) {

        if (myrand(g_rng) < plsOne) {
            // add +/- 1
            if (myrand(g_rng) < 0.5) {
                if (myrand(g_rng) < 0.5 || xcopy == 0) {
                    xcopy++;
                } else {
                    xcopy = xcopy - 1;
                }
            } else {
                if (myrand(g_rng) < 0.5 || ycopy == 0) {
                    ycopy++;
                } else {
                    ycopy = ycopy - 1;
                }
            }

        } else {
            if (myrand(g_rng) < 0.5) {
                if (myrand(g_rng) < 0.5 || xcopy == 1) {
                    xcopy = xcopy + 2;
                } else {
                    xcopy = xcopy - 2;
                }
            } else {
                if (myrand(g_rng) < 0.5 || ycopy == 1) {
                    ycopy = ycopy + 2;
                } else {
                    ycopy = ycopy - 2;
                }
            }
        }

    }
    _cnCalls[cell][segment] = std::make_pair(xcopy, ycopy);
}


//int SingleCell::gaussianDraw(int mean, double errorRate) {
//    //sample from a binomial distribution
//    boost::random::normal_distribution<double> distribution(mean, errorRate);
////    boost::math:: variate_generator<boost::mt19937, boost::math::normal_distribution<>> gen_norm(g_rng, distribution);
//    int draw = distribution(g_rng);
//    while (draw < 0) { //reject draws which are less than 0
//        draw = distribution(g_rng);
//    }
//    int copyNumber = std::round(draw);
//    return copyNumber;
//}

// _SCS_PREV[i][j] is the clone proportion of sample i for node j
// ecdf is the cumulative distribution for the clone proportions of the sample
// _nodes is the number of nodes/clones in the tree
void SingleCell::generateECDF(int i) {
    _ecdf.resize(1, std::vector<double>(_numClones, 0.0));
    double cdf = 0.0f;
    for (int j = 0; j < _numClones; ++j) {
        cdf += _SCS_PREV[j][i +
                            1]; //_SCS_PREV is a transpose of what ecdf should be ; the i+1 is because the first column of _SCS_PREV is the clone number
        _ecdf[0][j] = cdf;

    }

/*     for (int a = 0; a < _numClones; a++) {
         out << ecdf[0][a] << "\t";
     }  */

}

void SingleCell::initializeSCS() {

    int l = 0;
    //Getting total_mutations: 
    for (int i = 0; i < _nodeInformationRows; i++) {
        if (_NODE_INFORMATION[i][_mCol] > _numMutations) {
            _numMutations = _NODE_INFORMATION[i][_mCol];
        }
        l = i;
    }
    _numMutations++; //to account for the fact that there is a 0 mutation
    std::cout << "Total number of cells: " << _NUMSCS << std::endl;
    std::cout << "Total number of SNVs: " << _numMutations << std::endl;

    //ASSUMPTION: That there are all mutations consecutively (i.e., no missing mutation ids) or that it wouldn't matter if there is

    _varReads.resize(_NUMSCS, std::vector<int>(_numMutations));
    _refReads.resize(_NUMSCS, std::vector<int>(_numMutations));
    _totReads.resize(_NUMSCS, std::vector<int>(_numMutations));
    // _refAlleles.resize(_numClones, std::vector<std::vector<int>>(_numMutations, std::vector<int>(2, 0)));
    // _mutAlleles.resize(_numClones, std::vector<std::vector<int>>(_numMutations, std::vector<int>(2, 0))); 
    // _totalAlleles.resize(_numClones, std::vector<std::vector<int>>(_numMutations, std::vector<int>(2, 0))); 
    _segments.resize(_numClones, std::vector<int>(_numMutations));
    _segmentCopyNumbers.resize(_numClones, std::vector<std::vector<int>>(_numSegments, std::vector<int>(2, 0)));

    for (int i = 0; i < _NUMSCS; i++) {
        for (int j = 0; j < _numMutations; j++) {

            //_SCS[i][j] = 0;
            _varReads[i][j] = 0;
            _refReads[i][j] = 0;
            _totReads[i][j] = 0;

        }
       for(int ell=0; ell < _numSegments; ell++){
           _cnCalls[i][ell] = std::make_pair(0,0);
       }

    }
    std::cout << "done initializing read count matrices" << std::endl;
    for (int a = 0; a < _numClones; a++) {
        for (int b = 0; b < _numMutations; b++) {
            // _refAlleles[a][b][0] = 0;
            // _mutAlleles[a][b][0] = 0;
            // _totalAlleles[a][b][0] = 0;
            _segments[a][b] = 0;
            // _refAlleles[a][b][1] = 0;
            // _mutAlleles[a][b][1] = 0;
            // _totalAlleles[a][b][1] = 0;
        }
    }



    // for (int c = 0; c < _numClones; c++) {
    //     for (int d = 0; d < _numSegments; d++) {
    //         _segmentCopyNumbers[c][d][0] = 0;
    //         _segmentCopyNumbers[c][d][1] = 0;
    //     }
    // }



    for (int r = 0; r < _nodeInformationRows; r++) {
        int clone = _NODE_INFORMATION[r][_nodeCol];
        int seg = _NODE_INFORMATION[r][_segmentCol];
        int x = _NODE_INFORMATION[r][_xCol];
        int y = _NODE_INFORMATION[r][_yCol];



        _segmentCopyNumbers[clone][seg][0] = x;
        _segmentCopyNumbers[clone][seg][1] = y;
        // std::cout << clone << "," << seg << std::endl;

    }

    std::cout << "done initializing copy numbers" << std::endl;


}

//Generates variant and reference read counts for _NUMSCS single cells 
//sample is index of the biopsy sample and is used to make sure the correct clonal proporations are used
void SingleCell::generateCells(int sample, double cnaError, double plsOne) {
    std::cout << "generating single cell read counts" << std::endl;
    _cellLabels.resize(_NUMSCS);
    boost::random::uniform_01<double> myrand;
//    std::uniform_real_distribution<> myrand(0, 1);

    std::unordered_map<int, std::unordered_map<int, int>> mutationLookUp;
    for (int r = 0; r < _nodeInformationRows; r++) {
        int clone = _NODE_INFORMATION[r][_nodeCol];
        int mutation = _NODE_INFORMATION[r][_mCol];
        if (mutation != -1) {
            mutationLookUp[clone][mutation] = r;
        }
    }


    for (int i = 0; i < _NUMSCS; i++) { //i is the cell label

        //sample the clone id from the clonal proportions of the sample
        int clone = sampleSingleCells(sample);

        _cellLabels[i] = clone;


        // ASSUMPTION: that every clone will have every mutation

        for (int j = 0; j < _numMutations; j++) {

            //ASSUMPTION: that the node and mutation are present in the dataset
            int rowOfMutation = mutationLookUp[clone][j];
            if (!(rowOfMutation >= 0 && rowOfMutation < _nodeInformationRows)) {
                throw std::runtime_error("Error in generateCells: mutation lookup out of bounds");
            }
            assert(rowOfMutation >= 0 && rowOfMutation < _nodeInformationRows);

            int x = _NODE_INFORMATION[rowOfMutation][_xCol];
            int y = _NODE_INFORMATION[rowOfMutation][_yCol];

            int x_bar = _NODE_INFORMATION[rowOfMutation][_xbarCol];
            int y_bar = _NODE_INFORMATION[rowOfMutation][_ybarCol];





            //segment assignment doesn't change by clone, this can be a map
            _segments[clone][j] = _NODE_INFORMATION[rowOfMutation][_segmentCol];

            if (x + y > 0) {
                std::pair<int, int> exp = draw(x_bar + y_bar, x + y);
                _varReads[i][j] = exp.first;
                _totReads[i][j] = exp.second;
                _refReads[i][j] = _totReads[i][j] - _varReads[i][j];
            } else {
                _varReads[i][j] = 0;
                _totReads[i][j] = 0;
                _refReads[i][j] = _totReads[i][j] - _varReads[i][j];
            }




            //TO DO: Check how code handles this when there are no mutated alleles

        }
        //add copy number error

        for (int s = 0; s < _numSegments; s++) {
             addCNANoise(i, s, cnaError, plsOne );

        }
//            myFile << s << delim << i << delim << _segmentCopyNumbers[c][s][0] << delim
//                   << _segmentCopyNumbers[c][s][1] << std::endl;

        }

}


//sample a clone from the CDF of the clonal proportions for the specified sample
// FYI, there may be be some existing function to sample directly from the clonal props, but this works. 
// Feel free to replace if there is 
int SingleCell::sampleSingleCells(int sample) {

    float r;
    boost::random::uniform_01<> myrand;
//    std::uniform_real_distribution<> myrand(0, 1); //uniform distribution between 0 and 1

    r = myrand(g_rng);

    int index = 0;

    while (r > _ecdf[0][index] && index < _ecdf[0].size()) {
        index += 1;
    }

    if (index >= _ecdf[0].size()) {
        if (!(_ecdf[0][_ecdf[0].size() - 1] > .9999)) { //allows for .0001 in rounding error
            throw std::runtime_error("Error in sampleSingleCells: cumulative probability in ecdf is not 1");
        }
        assert (_ecdf[0][_ecdf[0].size() - 1] > .9999); //allows for .0001 in rounding error
        index = _ecdf[0].size() - 1;

    }

    return index;
}

// _alpha_fp is the per base sequencing error rate
// _READ_DEPTH is the specified coverage

std::pair<int, int> SingleCell::draw(int mut_alleles, int total_cn) {


    int treads;
    double cov = (_READ_DEPTH / 2) * (total_cn);


    boost::random::poisson_distribution<> readcounts(cov);
//    std::poisson_distribution<int> readcounts(cov);

    //draw total read counts 
    treads = readcounts(g_rng);
    int a = 0;


    //adjust success prob p for sequencing error
    double p = 1.0 * (mut_alleles) / (total_cn);
    p = p * (1 - _alpha_fp) + (1 - p) * _alpha_fp / 3;


    //draw variant read counts 
    int vreads = binomialdraw(p, treads);

    return (std::make_pair(vreads, treads));
}


int SingleCell::binomialdraw(double p, int n) {
    //sample from a binomial distribution
    boost::random::binomial_distribution<int> distribution(n, p);
    int draw = distribution(g_rng);

    return draw;
}

void SingleCell::printSCS(std::ostream &out, int sample) {


    std::string fname_sparse = _outFiles["SPARSE"] + ".p" + std::to_string(sample);
    std::string fname_cell = _outFiles["CELLS"] + ".p" + std::to_string(sample);
    std::string fname_cellAssignments = _outFiles["CELLASSIGNMENTS"] + ".p" + std::to_string(sample);




    //write a sparse dataframe of the variant and total read counts 
    if (_outFiles["SPARSE"].compare("") != 0) {
        std::cout << "Writing sparse file..." << std::endl;
        std::string delim = "\t";
        std::ofstream myFile(fname_sparse);
        myFile << "segment" << delim << "mutation" << delim << "cell" << delim << "varReads" << delim << "totReads"
               << std::endl;
        std::cout << "_numMutations: " << _numMutations <<std::endl;
        std::cout << "_NUMSCS: " << _NUMSCS <<std::endl;
        for (int j = 0; j < _numMutations; j++) {
            for (int i = 0; i < _NUMSCS; i++) {
                int seg = _segments[_cellLabels[i]][j];
//                 if(i == _NUMSCS -1 && j == _numMutations -1){
//                    std::cout << j << ":" << i << "seg: " << seg << " cellLabels[i] "  << _cellLabels[i] << ":"  << _varReads[i][j] << "," << _totReads[i][j] <<std::endl;
//
//                 }

                if (_totReads[i][j] > 0) {
                    myFile << seg << delim << j << delim << i << delim << _varReads[i][j] << delim << _totReads[i][j]
                           << std::endl;

                }
            }
        }
        myFile.close();
    }


    //write the cell labels 
    if (_outFiles["CELLS"].compare("") != 0) {
        std::string delim = ",";
        std::ofstream myFile(fname_cell);
        myFile << "segment" << delim << "cell" << delim << "copiesX" << delim << "copiesY" << std::endl;
        int prevSeg = -1;
        int seg = 0;
        for (int i = 0; i < _NUMSCS; i++) {
//            int c = _cellLabels[i];
            for (int s = 0; s < _numSegments; s++) {
                myFile << s << delim << i << delim << _cnCalls[i][s].first << delim
                       << _cnCalls[i][s].second << std::endl;

            }
        }
        myFile.close();
    }

    if (_outFiles["CELLASSIGNMENTS"].compare("") != 0) {
        std::string delim = ",";
        std::ofstream myFile(fname_cellAssignments);
        myFile << "Cell" << delim << "Cluster" << std::endl;
        for (int c = 0; c < _NUMSCS; c++) {
            myFile << c << delim << _cellLabels[c] << std::endl;
        }
        myFile.close();
    }

}

