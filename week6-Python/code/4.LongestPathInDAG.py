import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week6-python'
curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week6-python'
inputFile = './data/4.longest_path_in_DAG-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_245_7 (3).txt'
outputFile = './results/4.LongestPathInDAG.txt'


# set current directory
os.chdir(curDir)

# read input
with open(inputFile) as f:
    inputs = f.readlines()

source = int(inputs[0].strip())
sink = int(inputs[1].strip())
edgeInputs = [tuple([int(i) for i in re.split("->|:",s.strip())]) for s in inputs[2:len(inputs)]]

# create a graph using the edge inputs
g = nx.DiGraph()
g.add_weighted_edges_from(edgeInputs)

# set max length to negative infinite
NEGINF = float("-inf")
length = {n:NEGINF for n in nx.nodes(g)}
length[source] = 0
#print(length)

# set backtrack of each node, initially all empty
backtrack = {n:[] for n in nx.nodes(g)}

# create predecessor list of all the nodes
pred = {n:[] for n in nx.nodes(g)}
for u,v in nx.edges(g):
    pred[v] = pred[v] + [u]
pred[source] = []
print(pred)

## access all nodes in topological order and update length and backtrack
nodes = nx.topological_sort(g)
for n in nodes:
    print(n)
    if len(pred[n])>0:
        candidateLengths = [length[p]+g[p][n]['weight'] for p in pred[n]]
        length[n] = max(candidateLengths)
        if(length[n] > NEGINF):
            candidateIndex = candidateLengths.index(length[n])
            backtrack[n] = pred[n][candidateIndex]

## length of sink is the max length
longestLength = length[sink]


#print(nodes)

## find the path by backtracking from sink
longestPath = [sink]
n = sink
while n != source:
    print(backtrack[n])
    longestPath.append(backtrack[n])
    n = backtrack[n]
    
    
longestPath.reverse()
#print(longestPath)

# output
with open(outputFile, "w") as f:
    f.writelines(str(longestLength) + "\n")
    f.writelines("->".join([str(s) for s in longestPath]))

