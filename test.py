import clonelib as cl

# Example input values
# L = [(1, 1), (2, 1), (3,1)]
# L = [(1, 2), (2, 3), (3, 4)]
L = [(1,1), (1,0)]
root_x = 1  
root_y = 1  

# Call the function
result = cl.get_cna_trees(set(L), root_x, root_y)
print(f"Number of candidate CNA trees: {len(result)}")
for i in range(len(result)):
    trees = cl.get_genotype_trees(result[i])
    print(f"# of genotype trees: {len(trees)}")
    for j,tree in enumerate(trees):
        print(f"\nCNA tree: {i} Genotype tree {j}")
        print(tree)




  