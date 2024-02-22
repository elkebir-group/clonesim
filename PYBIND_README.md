# clonelib python interface
Python package to enumerate CNA and genotype trees from C++ code

Dependencies:
- python
- setuptools
- pybind11
- lemon 
- boost 

### Installation instructions:
modify the lemon and boost lib/include paths in `setup.py` as needed.
```
cd clonesim
python setup.py install
```

### Available Functions
Within the `clonelib` module there are two available functions: `get_cna_trees` and `get_genotype_trees`.

 - `get_cna_trees(L, root_x, root_y)` where `L` is a set of integer tuples representing CNA states (x,y) 
   while `root_x` and `root_y` is the root (root_x, root_y) of the CNA trees.
-

### Example usage
See `test.py` for example usage. 
```
import clonelib

# Example input values

L = [(1, 2), (2, 3), (3, 4)]
root_x = 1  
root_y = 1  

# Enumerate all possible CNA trees
result = clonelib.get_cna_trees(set(L), root_x, root_y)
print(f"Number of candidate CNA trees: {len(result)}")
for i in range(len(result)):
    trees = clonelib.get_genotype_trees(result[i])
    print(f"# of genotype trees: {len(trees)}")
    for j,tree in enumerate(trees):
        print(f"\CNA tree: {i} Genotype tree {j}")
        print(tree)


```
