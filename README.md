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

## Simulate

Simulate will generate the trees, proportions, and node information.

Inputs: 
- S : this is the input file of CNA trees.

Optional: 
- kk :  number of truncal segments
-  s : random generator seed
- purity: expected purity
- minProp: Minimum desired clone proportion
- n : Number of SNVs
- m: Number of samples
- k: Number of segments
- l: Number of mutation clusters
- dot: graphviz dot filename
- f : set if you would like to output the tree, node information, and proportion files
- output_file_dir: the directory you would like to put the output files in

 Example run: ./simulate -r -S /build/cnatrees.txt -purity .99 -minProp .05 -kk 2 -f -s 12 -l 7 -k 50 -n 5000 -m 1 -output_file_dir /build/output/intermediate 
 

## Python Interface
With the aid of pybind11, clonelib has a python interface to certain functionality within clonelib.
This functionality currently includes:  

    - enumeration of CNA trees give a set of allele-specific copy number states
    - enumeration of genotype trees given a CNA tree 

For more details on installation and usage, see [PYBIND_README.md](./PYBIND_README.md).


