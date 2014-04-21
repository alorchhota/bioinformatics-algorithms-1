import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
inputFile1 = './data/E-coli.txt'
inputFile2 = './data/Salmonella_enterica.fasta'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_289_5.txt'
outputFile = './results/4.SharedKmer.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile1) as f:
    inputs1 = f.readlines()

with open(inputFile2) as f:
    inputs2 = f.readlines()

k = 30
genome1 = inputs1[0].strip()
genome2 = ''.join([s.strip() for s in inputs2[1:]])

#k = 3
#genome1 = 'AAACTCATC'
#genome2 = 'TTTCAAATC'


del inputs1
del inputs2



print(len(genome1), len(genome2))

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

output = {}

n1 = len(genome1)
n2 = len(genome2)

print('getting g1kmers')
g1kmers = [genome1[i:i+k] for i in range(n1-k+1)]
genome1c = reverseComplement(genome1)
g1kmersc = [genome1c[n1-i-k:n1-i] for i in range(n1-k+1)]



## make reverse map of kmers to index
g2kmers = [genome2[i:i+k] for i in range(n2-k+1)]

print('generating inverted index')
print('init from g1kmers2')
g2kmerInvertedIndex = {kmer:[] for kmer in set(g2kmers)}

print('init from g1kmers')
g1kmerInvertedIndex = {kmer:[] for kmer in set(g1kmers)}
print('update')
g2kmerInvertedIndex.update(g1kmerInvertedIndex)
del g1kmerInvertedIndex

print('init from g1kmersc')
g1kmercInvertedIndex = {kmer:[] for kmer in set(g1kmersc)}
print('update')
g2kmerInvertedIndex.update(g1kmercInvertedIndex)
del g1kmercInvertedIndex

#for i in range(len(g2kmers)):
#    g2kmerInvertedIndex[g2kmers[i]] = []
print('finally, inverted index')
for i in range(len(g2kmers)):
    g2kmerInvertedIndex[g2kmers[i]].append(i)


print('finding points')
#print(g2kmerInvertedIndex)


#for i in range(n1-k+1):
    #kmer = genome1[i:i+k]
    #kmerc = genome1c[n1-i-k:n1-i]
    #ind = []
    #if kmer in g2kmersSet:
    #    ind = ind + g2kmerInvertedIndex[kmer]
    #if kmerc in g2kmersSet:
    #    ind = ind + g2kmerInvertedIndex[kmerc]

output = [set(g2kmerInvertedIndex[g1kmers[i]] + g2kmerInvertedIndex[g1kmersc[i]]) for i in range(n1-k+1)]


#g1kmers = [(genome1[i:i+k],genome1c[n1-i-k:n1-i]) for i in range(n1-k+1)]
#ind1 = [g2kmerInvertedIndex[g1kmers] for i in range(n1-k+1)]
        
        
with open(outputFile, 'w') as f:
    for key in range(len(output)):
        for x,y in zip([key]*len(output[key]),output[key]):
            f.writelines('(' + str(x) + ', ' + str(y) + ')\n')
    

# output
#with open(outputFile, 'w') as f:
#    f.writelines(str(d)+"\n")
