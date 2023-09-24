import networkx as nx
import pygraphviz as pgv

def read_tree(fname):
    tree = nx.DiGraph()
    with open(fname, "r+") as file:
        for edge in file:
    
            edge = edge.strip().split("\t")
            tree.add_edge(edge[0], edge[1])

    return tree 

def draw(tree, outfname):
    pgtree = pgv.AGraph( strict=False, directed=False)
    pgtree.add_edges_from(list(tree.edges))


    pgtree.layout("dot")
    pgtree.draw(outfname)

tpath = "/scratch/data/leah/pharming/simulation_study/clonesim/build/test"
fname = f"{tpath}/tree.txt"
tree = read_tree(fname)
draw(tree, f"{tpath}/tree.png")




