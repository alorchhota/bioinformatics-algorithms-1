import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
inputFile = './data/3.2BreakDistance-1.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_288_4.txt'
outputFile = './results/3.2BreakDistance.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome1 = inputs[0].strip()
genome2 = inputs[1].strip()

## parse input file and build genomes
startIndexes = [m.start() for m in re.finditer('[(]', genome1)]
endIndexes = [m.start() for m in re.finditer('[)]', genome1)]
g1 = [genome1[m[0]+1:m[1]] for m in zip(startIndexes, endIndexes)]
g1 = [[int(x) for x in re.split('\\s+', ch)] for ch in g1]

startIndexes = [m.start() for m in re.finditer('[(]', genome2)]
endIndexes = [m.start() for m in re.finditer('[)]', genome2)]
g2 = [genome2[m[0]+1:m[1]] for m in zip(startIndexes, endIndexes)]
g2 = [[int(x) for x in re.split('\\s+', ch)] for ch in g2]

## build 2-breakpoint graphs
bpg = nx.Graph()

## add connecting edges (except the edges)
# edge 1 has two ends: 1- to 1+ or (-1,+1)
# edge 2 has two ends: 2- to 2+ or (-2,+2) ... etc.
# chromosome (+1, +2, -3) has connections: (1+,2-); (2+,3+); (3-, 1-)

for ch in g1:
    n = len(ch)
    for i in range(n):
        bpg.add_edge(str(abs(ch[i]))+['-','+'][ch[i]>0],\
                     str(abs(ch[(i+1)%n]))+['+','-'][ch[(i+1)%n]>0])
    
for ch in g2:
    n = len(ch)
    for i in range(n):
        bpg.add_edge(str(abs(ch[i]))+['-','+'][ch[i]>0],\
                     str(abs(ch[(i+1)%n]))+['+','-'][ch[(i+1)%n]>0])


## number of connected components is the number of cycles
cc = nx.connected_components(bpg)
ncycles = len(cc)

## number of boxes is the total number of elements in genome
nboxes = sum([len(ch) for ch in g1])

## distance
d = nboxes - ncycles
print(d)

# output
with open(outputFile, 'w') as f:
    f.writelines(str(d)+"\n")
