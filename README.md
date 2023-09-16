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
./generatecnatrees > cnatrees.txt
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

## inputs
A list of flags to be set when calling simulate: 
(note : run ./simulate --help for most up-to date listing and default values)

Mandatory: 
- S : this is the input file of CNA trees.
- STree: Output filename for tree
- SNode: Output filename for csv with information about nodes, segments, and mutaitons
- SProportions: Output filename for csv with information about nodes, segments, and mutaitons
  

Optional: 
- kk :  number of truncal segments
- s : random generator seed
- purity: expected purity
- minProp: Minimum desired clone proportion
- n : Number of SNVs
- m: Number of samples
- k: Number of segments
- l: Number of mutation clusters
- e : Error rate for CNA data
- dot : Graphviz DOT output filename
- r : Remove unsampled nodes
- num_cells: The number of cells to simulate with the single cell generation
- read_depth: The read_depth for the single cell generation
- alpha_fp: The sequencing error for single cell generation
- out_dir: The output directory for single cell generation (Please note: single cell data will not output if this argument is missing)
- sc : Set to true to generate single cell data

## Example runs: 

To generate single cell data:

  With removing unsampled nodes: 
  
    ./simulate -n 50 -m 2 -purity .99 -SNode nodetest4.csv -STree treetest4.csv -SProportions ptest4.c -S cnatrees.txt  -alpha_fp .002 -r -sc -out_dir scdata
    
  Without removing unsampled nodes: 
  
    ./simulate -n 50 -m 2 -purity .99 -SNode nodetest4.csv -STree treetest4.csv -SProportions ptest4.c -S cnatrees.txt  -alpha_fp .002  -sc -out_dir testingtoseeifthisruns
    

Not generating single cell data:

  With removing unsampled nodes: 
  
    ./simulate -n 50 -m 2 -purity .99 -SNode nodetest4.csv -STree treetest4.csv -SProportions ptest4.c -S cnatrees.txt  -alpha_fp .002 -r
  
  Without removing unsampled nodes: 
  
    ./simulate -n 50 -m 2 -purity .99 -SNode nodetest4.csv -STree treetest4.csv -SProportions ptest4.c -S cnatrees.txt  -alpha_fp .002 



