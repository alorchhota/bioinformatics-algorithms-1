import os
import csv
import sys
import re
import importlib
import networkx as nx
import itertools

# settings
#curDir = 'E:/GitHub/bioinformatics-algorithms-1/week10'
curDir = 'D:/GitHub/bioinformatics-algorithms-1/week10'
inputFile = './data/2.inverse-bwt-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_297_4.txt'
outputFile = './results/2.inverse-bwt.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome = inputs[0].strip()

def ibwt(bwtgenome):
    """ Given BWT(genome), get genome. """
    firstRow = sorted([nt for nt in bwtgenome])
    print(''.join(firstRow))
    while 
    

genome = ibwt(genome)
print(genome)

## output
with open(outputFile, "w") as f:
    f.writelines(genome + '\n')

print('done.')
