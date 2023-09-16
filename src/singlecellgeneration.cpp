#include "singlecellgeneration.h"
#include <fstream>
#include <algorithm>
#include <string>
#include <lemon/connectivity.h>
#include <map>
#include <cmath>

SingleCell::SingleCell (int numSCS, float read_depth, float alpha_fp, std::string outdir, int numSegments, int numSamples, std::mt19937 &gen, double cna_error)
    : _NUMSCS(numSCS) //** HYPERPARAMETER
    , _READ_DEPTH(read_depth) //** HYPERPARAMETER
    , _gen(gen)
    , _numSegments(numSegments)
    , _numSamples(numSamples)
    , _cnaError(cna_error)
    , outdir(outdir)
    , _outFiles({ {"SPARSE", outdir + "/sparse"}, {"CELLS", outdir + "/cells"}, {"CELLASSIGNMENTS", outdir + "/cellAssignments"}})
    , _alpha_fp (alpha_fp)
    //, _seed
    , _varReads()
    , _refReads()
    , _totReads()
    , _cellLabels()
    , _SCS_PREV() //Note: the first column is the node ID 
    , _NODE_INFORMATION()
    , _ecdf() //cdf of clonal proportions
    , _mutAlleles()
    , _refAlleles ()
    , _totalAlleles()
    , _segmentCopyNumbers()
    , _segments()
    , _numClones (0)
    , _nodeInformationRows (0)
    , _nodeInformationColumns (0)
    , _numMutations (0)
    , _nodeCol (0)
    , _segmentCol (1)
    , _xCol (2)
    , _yCol (3)
    , _mCol (4)
    , _xbarCol (5)
    , _ybarCol (6)
    , _clusterIDCol (7)
    {

    }

void SingleCell::loadData(std::ostream&out, std::string&outputProportionFilename, std::string&outputNodeFilename) {
    int i = 0;
    int j = 0;
    std::ifstream file(outputProportionFilename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open node probability file");
    }
    std::string line; 
    getline(file, line); //discard header 
    std::string cell;
    float proportion = 0.0;
    while (std::getline(file, line)) {
        j = 0; 
        std::vector<float> row; 
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

    std::ifstream file2(outputNodeFilename);
    if (!file2.is_open()) {
        throw std::runtime_error("Coult not open node file");
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

    addNoiseCNA();

/*     for (int aa = 0; aa < _nodeInformationRows; aa++) {
        for (int bb = 0; bb < _nodeInformationColumns; bb++) {
            out << _NODE_INFORMATION[aa][bb] << "\t";
        }
        out << "\n";
    } 
 */
    

/* 
     for (int a = 0; a < i; a++) {
        for (int b = 0; b < j; b++) {
            out << _SCS_PREV[a][b] << "\t";
        }
        out << "\n"; */
    //}  
    
    /* for (int aa = 0; aa < _nodeInformationRows; aa++) {
        for (int bb = 0; bb < _nodeInformationColumns; bb++) {
            out << _NODE_INFORMATION[aa][bb] << "\t";
        }
        out << "\n";
    }  */
    return;
}

void SingleCell::addNoiseCNA() {
    for (int r = 0; r < _nodeInformationRows; r++) {
        int xcopy = _NODE_INFORMATION[r][_xCol];
        int ycopy = _NODE_INFORMATION[r][_yCol];
        int newX = gaussianDraw(xcopy, _cnaError);
        int newY = gaussianDraw(ycopy, _cnaError);
        _NODE_INFORMATION[r][_xCol] = newX;
        _NODE_INFORMATION[r][_yCol] = newY;
    }
}

int SingleCell::gaussianDraw(int mean, double errorRate)
{
    //sample from a binomial distribution
    double stdDev = errorRate*mean;
    std::normal_distribution<double> distribution(mean, stdDev);
    int draw = distribution(_gen);
    if (draw < 0) {
        draw = 0;
    }
    int copyNumber = std::round(draw);
    return copyNumber;
}

// _SCS_PREV[i][j] is the clone proportion of sample i for node j
// ecdf is the cumulative distribution for the clone proportions of the sample
// _nodes is the number of nodes/clones in the tree
void SingleCell::generateECDF(std::ostream&out, int i) 
{
    _ecdf.resize(1, std::vector<float>(_numClones, 0.0));
    float cdf = 0.0f;
    for (int j = 0; j < _numClones; ++j)
    {
        cdf += _SCS_PREV[j][i+1]; //_SCS_PREV is a transpose of what ecdf should be ; the i+1 is because the first column of _SCS_PREV is the clone number
        _ecdf[0][j] = cdf;
 
    }

/*     for (int a = 0; a < _numClones; a++) {
         out << ecdf[0][a] << "\t";
     }  */

}

void SingleCell::initializeSCS(std::ostream&out)
{

    //Getting total_mutations: 
    for (int i = 0; i < _nodeInformationRows; i++) {
        if (_NODE_INFORMATION[i][_mCol] > _numMutations) {
            _numMutations = _NODE_INFORMATION[i][_mCol];
        }
    }
    _numMutations ++; //to account for the fact that there is a 0 mutation
    
    //ASSUMPTION: That there are all mutations consecutively (i.e., no missing mutation ids) or that it wouldn't matter if there is

    _varReads.resize(_NUMSCS, std::vector<int>(_numMutations));;
    _refReads.resize(_NUMSCS, std::vector<int>(_numMutations));;
    _totReads.resize(_NUMSCS, std::vector<int>(_numMutations));;
    _refAlleles.resize(_numClones, std::vector<std::vector<int>>(_numMutations, std::vector<int>(2, 0)));; 
    _mutAlleles.resize(_numClones, std::vector<std::vector<int>>(_numMutations, std::vector<int>(2, 0)));; 
    _totalAlleles.resize(_numClones, std::vector<std::vector<int>>(_numMutations, std::vector<int>(2, 0)));; 
    _segments.resize(_numClones, std::vector<int>(_numMutations));
    _segmentCopyNumbers.resize(_numClones, std::vector<std::vector<int>>(_numSegments, std::vector<int>(2, 0)));;
    for (int i = 0; i < _NUMSCS; i++) 
    {
        for (int j = 0; j < _numMutations; j++) 
        {

            //_SCS[i][j] = 0;
            _varReads[i][j] = 0;
            _refReads[i][j] = 0;
            _totReads[i][j] = 0;
    
        }
    }

    for (int a = 0; a < _numClones; a++) {
        for (int b = 0; b < _numMutations; b++) {
            _refAlleles[a][b][0] = 0;
            _mutAlleles[a][b][0] = 0;
            _totalAlleles[a][b][0] = 0;
            _segments[a][b] = 0;
            _refAlleles[a][b][1] = 0;
            _mutAlleles[a][b][1] = 0;
            _totalAlleles[a][b][1] = 0;
        }
    }

    for (int c = 0; c < _numClones; c++) {
        for (int d = 0; d < _numSegments; d++) {
            _segmentCopyNumbers[c][d][0] = 0;
            _segmentCopyNumbers[c][d][1] = 0;
        }
    }

    for (int r = 0; r < _nodeInformationRows; r++) {
        int clone = _NODE_INFORMATION[r][_nodeCol];
        int seg = _NODE_INFORMATION[r][_segmentCol];
        int x = _NODE_INFORMATION[r][_xCol];
        int y = _NODE_INFORMATION[r][_yCol];
        _segmentCopyNumbers[clone][seg][0] = x;
        _segmentCopyNumbers[clone][seg][1] = y;
    }
}

//Generates variant and reference read counts for _NUMSCS single cells 
//sample is index of the biopsy sample and is used to make sure the correct clonal proporations are used
void SingleCell::generateCells(std::ostream&out, int sample){
    _cellLabels.resize(_NUMSCS);
    std::uniform_real_distribution<> myrand(0, 1);
    for(int i=0; i < _NUMSCS; i++){ //i is the cell label

        //sample the clone id from the clonal proportions of the sample
        int clone = sampleSingleCells(std::cout, sample);
        
        _cellLabels[i] = clone;


            // ASSUMPTION: that every clone will have every mutation
            for (int j = 0; j < _numMutations; j++)
            {
                

                //Looping through node information to look for the clone and mutation: 
                //ASSUMPTION: that the node and mutation are present in the dataset; if not, this will attempt to index at -1
                int rowOfMutation = -1; 
                for (int r = 0; r < _nodeInformationRows; r++) {
                    if (_NODE_INFORMATION[r][_nodeCol] == clone && _NODE_INFORMATION[r][_mCol] == j) {
                        rowOfMutation = r;
                        break;
                    }
                }

                _totalAlleles[clone][j][0] = _NODE_INFORMATION[rowOfMutation][_xCol];
                _totalAlleles[clone][j][1] = _NODE_INFORMATION[rowOfMutation][_yCol];

                _mutAlleles[clone][j][0] = _NODE_INFORMATION[rowOfMutation][_xbarCol];
                _mutAlleles[clone][j][1] = _NODE_INFORMATION[rowOfMutation][_ybarCol];

                _refAlleles[clone][j][0] = _totalAlleles[clone][j][0] - _mutAlleles[clone][j][0];
                _refAlleles[clone][j][1] = _totalAlleles[clone][j][1] - _mutAlleles[clone][j][1];

                _segments[clone][j] = _NODE_INFORMATION[rowOfMutation][_segmentCol];

                std::pair<int, int> exp = draw(out, _mutAlleles[clone][j][0] +  _mutAlleles[clone][j][1], _refAlleles[clone][j][0] + _refAlleles[clone][j][1]);


            
                _varReads[i][j] = exp.first;
                _totReads[i][j] = exp.second; 
                _refReads[i][j] = _totReads[i][j] - _varReads[i][j];


                //TO DO: Check how code handles this when there are no mutated alleles

            }
    }
}


//sample a clone from the CDF of the clonal proportions for the specified sample
// FYI, there may be be some existing function to sample directly from the clonal props, but this works. 
// Feel free to replace if there is 
int SingleCell::sampleSingleCells(std::ostream&out, int sample)
{
    float r;
    std::uniform_real_distribution<> myrand(0, 1); //uniform distribution between 0 and 1
    r = myrand(_gen);

    int index = 0;
    float total = 0;

    while (r > _ecdf[0][index])
    {
        index += 1;
    }

    return index;
}

// _alpha_fp is the per base sequencing error rate
// _READ_DEPTH is the specified coverage

std::pair<int, int> SingleCell::draw(std::ostream&out, int mut_alleles, int ref_alleles)
{

 
    int treads;
    float cov = (_READ_DEPTH/2) * (mut_alleles + ref_alleles);


    std::poisson_distribution<int> readcounts(cov);
    
    //draw total read counts 
    treads = readcounts(_gen);
    int a = 0;


    //adjust success prob p for sequencing error
    float p = 1.0 * (mut_alleles) / (mut_alleles + ref_alleles);
    p = p * (1 - _alpha_fp) + (1 - p) * _alpha_fp/3;


    //draw variant read counts 
    int vreads = binomialdraw(p, treads);

    return (std::make_pair(vreads, treads));
}


int SingleCell::binomialdraw(float p, int n)
{
    //sample from a binomial distribution
    std::binomial_distribution<int> distribution(n, p);
    int draw = distribution(_gen);

    return draw;
}

void SingleCell::printSCS(std::ostream&out, int sample)
{


    std::string fname_sparse = _outFiles["SPARSE"] + ".p" + std::to_string(sample);
    std::string fname_cell = _outFiles["CELLS"] + ".p" + std::to_string(sample);
    std::string fname_cellAssignments = _outFiles["CELLASSIGNMENTS"] + ".p" + std::to_string(sample);



 
    //write a sparse dataframe of the variant and total read counts 
    if (_outFiles["SPARSE"].compare("") != 0)
    {
        std::string delim = "\t";
        std::ofstream myFile(fname_sparse);
        myFile << "segment" << delim << "mutation" << delim << "cell" << delim << "varReads" << delim << "totReads" << std::endl; 

        for (int j = 0; j < _numMutations; j++)
        {
            for (int i = 0; i < _NUMSCS; i++)
            {
                int seg = _segments[_cellLabels[i]][j];
    
                if (_totReads[i][j] > 0)
                {
                    myFile<< seg << delim << j << delim << i << delim  << _varReads[i][j] << delim << _totReads[i][j] << delim << std::endl;
                
                }
            }
        }
        myFile.close();
    }


    //write the cell labels 
    if (_outFiles["CELLS"].compare("") != 0)
    {
        std::string delim = ",";
        std::ofstream myFile(fname_cell);
        myFile << "segment" << delim << "cell" << delim << "copiesX" << delim << "copiesY" << std::endl;
        int prevSeg = -1;
        int seg = 0;
        for (int i = 0; i < _NUMSCS; i++)
        {
            int c = _cellLabels[i];
            for (int s = 0; s < _numSegments; s++) {
                myFile << s << delim << i << delim << _segmentCopyNumbers[c][s][0] << delim << _segmentCopyNumbers[c][s][1] << std::endl;
                
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

