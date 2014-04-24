import os
import csv
import sys
import re
import importlib
import networkx as nx
import itertools

# settings
curDir = 'E:/GitHub/bioinformatics-algorithms-1/week10'
#curDir = 'D:/GitHub/bioinformatics-algorithms-1/week10'
inputFile = './data/2.inverse-bwt-2.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_299_10.txt'
outputFile = './results/2.inverse-bwt.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome = inputs[0].strip()

def ibwt(bwtgenome):
    """ Given BWT(genome), get genome. """

    genomeLen = len(bwtgenome)
    
    # construct the first col
    firstCol = sorted([nt for nt in bwtgenome])

    # find unique letters in the genome
    alphabets = sorted(list(set(firstCol)))

    # find the index of first occurence of each letter in the first col
    start, firstColIndex = 0, {}
    for letter in alphabets:
        firstColIndex[letter] = start = firstCol.index(letter, start)
        
    # count the occurence in the last column
    freq = {letter:0 for letter in alphabets}
    lastColCount = [0]*genomeLen
    for i in range(genomeLen):
        letter = bwtgenome[i]
        lastColCount[i] = freq[letter]
        freq[letter] += 1
        
    # iterate using first-last property
    fidx = 0
    genome = ['$']
    while True:
        letter = bwtgenome[fidx]
        if letter == '$':
            break
        genome.append(letter)
        fidx = firstColIndex[letter] + lastColCount[fidx]

    # reversed genome is the original genome
    genome.reverse()
    return ''.join(genome)

genome = ibwt(genome)
#print(genome)

## output
with open(outputFile, "w") as f:
    f.writelines(genome + '\n')

print('done.')
