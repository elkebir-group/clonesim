//definitions

const int _NUMSCS;   //number of single cells to sample
const float _READ_DEPTH;   //sequencing coverage

std::mt19937 gen;  //rng generator, make sure you are specifying a seed 
 // gen(_seed)

//these are for a single sample, add a dimension if you want a sample index
std::vector<std::vector<int>> _varReads;
std::vector<std::vector<int>> _refReads;
std::vector<int> _cellLabels;   //indexed by cell, stores the clone id of the sampled cell 

//index by sample and clone/node to access the clonal proportions
std::vector<std::vector<float>>_SCS_PREV;
std::vector<std::vector<float>> ecdf;  //stores the CDF of the clonal proportions



std::vector<std::vector<int>> _mutAlleles;  //mutated copies
std::vector<std::vector<int>> _refAlleles; //reference copies


// First generateECDF is called to convert the clonal proportions to one that can be sampled from
// FYI, there may be be some existing function to sample directly from the clonal props, but this works. Feel free to replace

// Second initializeSCS just clears _varReads, _refReads, and _cellLabels, you don't need this if you add an inded for sample

// Third: generateCells generates the read counts and ground truth cell samples
// Fourth: printSCS writes the files to disk



   std::cout << "Generating Single Cell Data..." << std::endl;


  

        for (int i = 0; i < _samples; i++)
        {
            generateECDF(i);
            initializeSCS();

            generateCells(i);

            std::cout << "Writing Single Cell Data to Disk..." << std::endl;
            printSCS(i);
        }



void Simulate::initializeSCS()
{

    for (int i = 0; i < _NUMSCS; i++)
    {
        for (int j = 0; j < _total_mutations; j++)
        {

            _SCS[i][j] = 0;
            _varReads[i][j] = 0;
            _refReads[i][j] = 0;
    
        }
    }
}

// _SCS_PREV[i][j] is the clone proportion of sample i for node j
// ecdf is the cumulative distribution for the clone proportions 
// _nodes is the number of nodes/clones in the tree
void Simulate::generateECDF(int i)
{

    float cdf = 0.0f;
    for (int j = 0; j < _nodes; ++j)
    {
        cdf += _SCS_PREV[i][j];
        ecdf[i][j] = cdf;
 
    }
}



//Generates variant and reference read counts for _NUMSCS single cells 
//sample is index of the biopsy sample and is used to make sure the correct clonal proporations are used
void Simulate::generateCells(int sample){

    std::uniform_real_distribution<> myrand(0, 1);
    for(int i=0; i < _NUMSCS; i++){


        //sample the clone id from the clonal proportions of the sample
        int clone = sampleSingleCells(sample);
        
        _cellLabels[i] = clone;
       
  
    
            for (int j = 0; j < _total_mutations; j++)
            {
                

                //replace the below code with code to access the genotype (x,y, x_bar, y_bar) for SNV j of the sampled clone
                // int mut_bin = _snvBin[j];

                //_totAlleles = x + y = w
                // _mutAlleles = max(x_bar, y_bar)
                // _refAlleles =  _totaAlleles - _mutAlleles

                // int mut_alleles = _mutAlleles[clone][j];
                // int ref_alleles = _refAlleles[clone][j];   
                // end to modify

                std::pair<int, int> exp = draw(mut_alleles, ref_alleles);

            
                _varReads[i][j] = exp.first;
                
        
                
                _totReads[i][j] = exp.second 

            }



    }

}

//sample a clone from the CDF of the clonal proportions for the specified sample
// FYI, there may be be some existing function to sample directly from the clonal props, but this works. 
// Feel free to replace if there is 
int Simulate::sampleSingleCells(int sample)
{

    std::uniform_real_distribution<> myrand(0, 1); //uniform distribution between 0 and 1
    float r = myrand(gen);

    int index = 0;
    float total = 0;

    while (r > ecdf[sample][index])
    {
        index += 1;
    }

    return index;
}



// _alpha_fp is the per base sequencing error rate, make default 0.001
// _READ_DEPTH is the specified coverage

std::pair<int, int> Simulate::draw(int mut_alleles, int ref_alleles)
{

 
    int treads;

    float cov = (_READ_DEPTH/2) * (mut_alleles + ref_alleles);
    std::poisson_distribution<int> readcounts(cov);
    
    //draw total read counts 
    treads = readcounts(gen);

    //adjust success prob p for sequencing error
    float p = 1.0 * (mut_alleles) / (mut_alleles + ref_alleles);
    p = p * (1 - _alpha_fp) + (1 - p) * _alpha_fp/3;


    //draw variant read counts 
    int vreads = binomialdraw(p, treads, gen);

    return (std::make_pair(vreads, treads));
}


static int binomialdraw(float p, int n, std::mt19937 gen)
{
    //sample from a binomial distribution
    std::binomial_distribution<int> distribution(n, p);
    int draw = distribution(gen);

    return draw;
}


void Simulate::printSCS(int sample)
{

   
    std::string fname_alt = _outFiles["AD"];
    std::string fname_sparse = _outFiles["SPARSE"];
    std::string fname_tot = _outFiles["DP"];



    if (_samples > 1)
    {
        std::string fname_csv = _outFiles["SCS"] + ".p" + std::to_string(sample);
        std::string fname_alt = _outFiles["AD"] + ".p" + std::to_string(sample);
        std::string fname_tot = _outFiles["DP"] + ".p" + std::to_string(sample);

    }





 
    //write a sparse dataframe of the variant and total read counts 
    if (_outFiles["SPARSE"].compare("") != 0)
    {
        std::string delim = "\t";
        std::ofstream myFile(fname_sparse);

        for (int j = 0; j < _total_mutations; j++)
        {
            for (int i = 0; i < _NUMSCS; i++)
            {
                int total = _varReads[i][j] + _refReads[i][j];
                int copies = _refAlleles[_cellLabels[i]][j] + _mutAlleles[_cellLabels[i]][j];
            
                if (total > 0)
                {
                    std::pair<int, int> res = lookup_chrom_arm_by_bin(_snvBin[j]);
                    myFile <<  res.first << delim << j << delim << i << delim <<"G" << delim  << _varReads[i][j] << delim << total <<delim << copies << std::endl;
                }
            }
        }
        myFile.close();
    }


    //write the cell labels 
    if (_outFiles["CELLS"].compare("") != 0)
    {
        std::string delim = ",";
        std::ofstream myFile(_outFiles["CELLS"]);
        myFile << "cell" << delim << "cluster" << std::endl;
        for (int i = 0; i < _NUMSCS; i++)
        {

            myFile << i << delim << _cellLabels[i] << std::endl;
        }
        myFile.close();
    }



  
    //write wide matrix of variation read counts 
    if (_outFiles["AD"].compare(""))
    {
        write_csv(fname_alt, _cellIDs, _colnames, _varReads, ",");
    }


    //write wide matrix of total read counts 
    if (_outFiles["DP"].compare("") != 0)
    {
        write_csv(fname_tot, _cellIDs, _colnames, _refReads, ",");
    }

}



//https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/
void Simulate::write_csv(std::string filename, std::vector<std::string> cell_ids,
                         std::vector<std::string> colnames,
                         std::vector<std::vector<int> > dataset, std::string delim)
{

    // Create an output filestream object
    std::ofstream myFile(filename);
    // Send column names to the stream
    for (int j = 0; j < colnames.size(); ++j)
    {

        myFile << colnames[j];
        if (j != colnames.size() - 1)
            myFile << delim; // No comma at end of line
    }
    if (colnames.size() > 0)
    {
        myFile << "\n";
    }

    // Send data to the stream
    for (int i = 0; i < dataset.size(); ++i)
    {
        myFile << cell_ids[i] << delim;
        for (int j = 0; j < dataset[i].size(); ++j)
        {
            myFile << dataset[i][j];
            if (j != dataset[i].size() - 1)
                myFile << delim; // No comma at end of line
        }
        myFile << "\n";
    }

    // Close the file
    myFile.close();
}