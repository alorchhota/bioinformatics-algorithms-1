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
inputFile = './data/6.MultiplePatternMatching-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_301_7.txt'
outputFile = './results/6.MultiplePatternMatching.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

bwtgenome = inputs[0].strip()
patterns = re.split('\\s+', inputs[1].strip())

#print(genome)
#print(patterns)
#bwtgenome = 'smnpbnnaaaaa$a'
genomeLen = len(bwtgenome)

# construct the first col
firstCol = sorted([nt for nt in bwtgenome])

# find unique letters in the genome
alphabets = sorted(list(set(firstCol)))

# find the index of first occurence of each letter in the first col
start, firstOccurence = 0, {}
for letter in alphabets:
    firstOccurence[letter] = start = firstCol.index(letter, start)
    
# count the occurence in the last column and build the count matrix
#freq = {letter:0 for letter in alphabets}
count = {letter:[0]*(genomeLen+1) for letter in alphabets}
for i in range(genomeLen):
    for letter in alphabets:
        count[letter][i+1] = count[letter][i]
    count[bwtgenome[i]][i+1] += 1
    
#print(count)

def bwtMatch(pattern):
    top = 0
    bottom = genomeLen - 1

    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern[-1]
            pattern = pattern[0:-1]
            if count[symbol][bottom+1] - count[symbol][top] > 0:
                top = firstOccurence[symbol] + count[symbol][top]
                bottom = firstOccurence[symbol] + count[symbol][bottom+1] - 1
            else:
                return 0
        else:
            return bottom - top + 1

matchCount = [bwtMatch(p) for p in patterns]
#print(matchCount)

## output
with open(outputFile, "w") as f:
    f.writelines(' '.join([str(n) for n in matchCount]) + '\n')

print('done.')
