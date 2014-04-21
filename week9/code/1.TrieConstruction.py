import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
curDir = 'E:/GitHub/bioinformatics-algorithms-1/week9'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
inputFile = './data/1.TrieConstruction-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_294_3.txt'
outputFile = './results/1.TrieConstruction.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

patterns = [s.strip() for s in inputs]
#print(patterns)

## trie graph
g = nx.DiGraph()
g.add_node(1)     # add root with id 1

#g.add_edges_from([(1,2,{'label':'T'}),(1,3,{'label':'A'}),(1,4,{'label':'C'}),(2,5,{'label':'A'})])

neighborsWithLabel = lambda node, label: [e[1] for e in g.edges_iter(node, data=True) if e[2]['label']==label]
getNewNode = lambda : len(g.nodes())+1

#print(neighborsWithLabel(1,'A'))
#print(getNewNode())

for p in patterns:
    # traverse as long as pattern matches
    curNode = 1
    for i in range(len(p)):
        nextNode = neighborsWithLabel(curNode, p[i])
        if len(nextNode) == 0:
            break
        curNode = nextNode[0]

    # add edges for the remaining unmatched pattern
    if len(nextNode) == 0:
        for j in range(i, len(p)):
            newNode = g.number_of_nodes()+1
            g.add_edge(curNode, newNode, {'label':p[j]})
            curNode = newNode
    
adj = [[e[0], e[1], e[2]['label']] for e in g.edges_iter(data=True)]


## output
with open(outputFile, "w") as f:
    f.writelines('\n'.join([' '.join([str(s) for s in row]) for row in adj]))
