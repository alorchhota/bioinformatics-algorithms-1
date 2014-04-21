import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
inputFile = './data/4.SharedKmer-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_289_5.txt'
outputFile = './results/4.SharedKmer.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

k = int(inputs[0].strip())
genome1 = inputs[1].strip()
genome2 = inputs[2].strip()

def reverseComplement(genome):
    rcmap = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rc = [rcmap[nt] for nt in genome]
    rc.reverse()
    return ''.join(rc)

def complement(genome):
    rcmap = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rc = [rcmap[nt] for nt in genome]
    return ''.join(rc)


print(len(genome1))
genome1c = reverseComplement(genome1)

output = {}

def writeOutput(f):
    global output
    s = ''
    #print(output.keys())
    for key in output.keys():
        for x,y in zip([key]*len(output[key]),output[key]):
            s += '(' + str(x) + ',' + str(y) + ')\n'
    #print(s)
    f.writelines(s)
    output.clear()


n1 = len(genome1)
n2 = len(genome2)

## make reverse map of kmers to index
g2kmers = [genome2[i:i+k] for i in range(n2-k+1)]
g2kmersSet = set(g2kmers)
g2kmerInvertedIndex = {kmer:[] for kmer in g2kmersSet}
for i in range(len(g2kmers)):
    g2kmerInvertedIndex[g2kmers[i]].append(i)

#print(g2kmerInvertedIndex)

for i in range(n1-k+1):
    kmer = genome1[i:i+k]
    kmerc = genome1c[n1-i-k:n1-i]
    ind = []
    if kmer in g2kmersSet:
        ind = ind + g2kmerInvertedIndex[kmer]
    if kmerc in g2kmersSet:
        ind = ind + g2kmerInvertedIndex[kmerc]
    output[i] = set(ind)
        
        
with open(outputFile, 'w') as f:
    for key in output.keys():
        for x,y in zip([key]*len(output[key]),output[key]):
            f.writelines('(' + str(x) + ', ' + str(y) + ')\n')
    

# output
#with open(outputFile, 'w') as f:
#    f.writelines(str(d)+"\n")
