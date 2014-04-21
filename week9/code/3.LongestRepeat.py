import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
curDir = 'E:/GitHub/bioinformatics-algorithms-1/week9'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
inputFile = './data/3.LongestRepeat-2.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_294_7.txt'
outputFile = './results/3.LongestRepeat.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome = inputs[0].strip() + '$'
#print(genome)
#print(patterns)

## build trie graph
g = nx.DiGraph()
g.add_node(1)     # add root with id 1

# two required function
neighborsWithLabel = lambda node, label: [e[1] for e in g.edges_iter(node, data=True) if e[2]['label']==label]
getNewNode = lambda : len(g.nodes())+1

longestRepeat = '';
maxLength = 0;

for idx in range(len(genome)):
    print(idx)
    p = genome[idx:]
    
    # traverse as long as pattern matches
    curNode = 1
    for i in range(len(p)):
        nextNode = neighborsWithLabel(curNode, p[i])
        if len(nextNode) == 0:
            break
        curNode = nextNode[0]

    length = i
    if length > maxLength:
        maxLength = length
        longestRepeat = genome[idx:idx+i]

    # add edges for the remaining unmatched pattern
    if len(nextNode) == 0:
        for j in range(i, len(p)):
            newNode = g.number_of_nodes()+1
            g.add_edge(curNode, newNode, {'label':p[j]})
            curNode = newNode

        # add the index at last
        g.node[curNode]['start'] = idx

print(longestRepeat)

#### find match
##def matchGenomeAt(idx):
##    match = False
##    curNode = 1
##    n = len(genome)
##    while(True):
##        # does not match, if idx >= len(genome)
##        if idx >= n:
##            break
##        
##        # does not match is no edge is found
##        nextNode = neighborsWithLabel(curNode, genome[idx])
##        if len(nextNode) == 0:
##            break
##
##        # move to next node in the tree
##        curNode = nextNode[0]
##        idx = idx + 1
##
##        # match if 'match' attribute of node is True
##        if 'match' in g.node[curNode].keys() and g.node[curNode]['match']==True:
##            match = True
##            break
##        
##    return match
##    
##matchIndexes = [i for i in range(len(genome)) if matchGenomeAt(i)==True]
###print(matchIndexes)    
##
##
#### output
##with open(outputFile, "w") as f:
##    f.writelines(' '.join([str(m) for m in matchIndexes]))
##
print('done.')
