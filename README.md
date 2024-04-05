# clonesim

## compilation

```
mkdir build
cd build
cmake ../ -DLIBLEMON_ROOT=/scratch/software/lemon/  -DBOOST_ROOT=/scratch/software/boost_1_74_0/
make
```

## usage

```
./generatecnatrees -minTotalCN 1 > cnatrees.txt
./simulate -S cnatrees.txt -dot T.dot > T.txt
dot -Tpdf T.dot > T.pdf
```

T.txt contains a text representantation of the simulated clonal tree. This can be visualized using:

```
./visualizephylo T.txt | dot -Tpdf > T.pdf
```

To remove unsampled nodes, use `-r` in `simulate`:


```
./simulate -r -n 50 -S cnatrees.txt -dot T.dot > T.txt
dot -Tpdf T.dot > T.pdf
```

## Simulate

Simulate will generate the trees, proportions, and node information.

```
Usage:
  ./simulate [--help|-h|-help] -S str [-dirich_param num] [-dot str] [-f]
     [-k int] [-kk int] [-l int] [-m int] [-minProp num] [-n int]
     [-num_tries int] [-output_file_dir str] [-purity num] [-r]
     [-restrictLoss] [-s int] [-threshold num] [-uniform]
Where:
  --help|-h|-help
     Print a short help message
  -S str
     Input CNA tree file
  -dirich_param num
     symmetric concentration parameter for Dirichlet distribution (default 2)
  -dot str
     Graphviz DOT output filename (default: '', no output)
  -f
     Whether to output files
  -k int
     Number of segments (default: 10)
  -kk int
     Number of truncal segments (default: 5)
  -l int
     Number of mutation clusters (default: 5)
  -m int
     Number of samples (default: 2)
  -minProp num
     Minimum desired clone proportion (default: 0.05)
  -n int
     Number of SNVs (default: 1000)
  -num_tries int
     The number of tries for sampling mutation rejection sampling (default 1000)
  -output_file_dir str
     The directory for where to write output files
  -purity num
     Expected purity (default: 0.8)
  -r
     Remove unsampled nodes
  -restrictLoss
     Whether to restrict copy number loss (default false)
  -s int
     Random number generator seed (default: 0)
  -threshold num
     Minimum threshold for SNV proportions (default: 0.05)
  -uniform
     use uniform distribution for mutation assignments
```



 Example run:
 ```
  ./simulate -r -S /build/cnatrees.txt -purity .99 -minProp .05 -kk 2 -f -s 12 -l 7 -k 50 -n 5000 -m 1 -output_file_dir /build/output/intermediate 
```

## Generate Single Cells

Generate single cells will generate single cell data for the simulation instance

Usage: ./generatesinglecells [-num_cells] [-read_depth] [-alpha_fp] [-out_dir] -in_dir [k] [-m] [-e] [-add] [-sub] 
Where: 
    -num_cells The number of cells to simulate with the single cell generation (default: 1000)
    -read_depth The read_depth for the single cell generation (default: .05)
    -alpha_fp The sequencing error for single cell generation (default .001)
    -out_dir The output directory for single cell generation (default: results)
    -in_dir The input directory of files for single cell generation (mandatory)
    -k Number of segments (default: 10))
     -m Number of samples (default: 2)
     -e The error rate for CNA data (default: 0)
      -add The proportion of copy number errors adding an allele (default .3
      -sub The proportion of copy number errors subtracting an allele (default .3)

      (Please note: total copy number neutral CNAs will have a probability of 1 - (add + sub) )

Example run: ./generatesinglecells -num_cells 500 -read_depth 0.05 -alpha_fp 0.001 -k 50 -m 1 -e .3 -add .35 -sub .25 -out_dir /Users/annahart/CLionProjects/clonesim3/build/output/testing -in_dir /Users/annahart/CLionProjects/clonesim3/build/output/intermediate


## Python Interface
With the aid of pybind11, clonelib has a python interface to certain functionality within clonelib.
This functionality currently includes:  

- enumeration of CNA trees give a set of allele-specific copy number states
- enumeration of genotype trees given a CNA tree 

For more details on installation and usage, see [PYBIND_README.md](./PYBIND_README.md).


