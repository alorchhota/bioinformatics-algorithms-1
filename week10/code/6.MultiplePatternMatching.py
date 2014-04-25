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
#inputFile = 'C:/Users/Ashis/Downloads/dataset_303_4.txt'
outputFile = './results/6.MultiplePatternMatching.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome = inputs[0].strip() + '$'
patterns = [p.strip() for p in inputs[1:]]


print(len(genome))
#print(len(patterns[0]))
#genome = genome[1:9000]
#patterns = patterns[0:5]

#print(genome)
#print(patterns)

def suffixArray(genome):
    suffix  = lambda idx: genome[idx:]
    
    def merge_sort(m):
        if len(m) <= 1:
            return m
     
        middle = int(len(m) / 2)
        left = m[:middle]
        right = m[middle:]
     
        left = merge_sort(left)
        right = merge_sort(right)
        return list(merge(left, right))


    def merge(left, right):
        result = []
        left_idx, right_idx = 0, 0
        while left_idx < len(left) and right_idx < len(right):
            # change the direction of this comparison to change the direction of the sort
            if suffix(left[left_idx]) <= suffix(right[right_idx]):
                result.append(left[left_idx])
                left_idx += 1
            else:
                result.append(right[right_idx])
                right_idx += 1
     
        if left:
            result.extend(left[left_idx:])
        if right:
            result.extend(right[right_idx:])
        return result

    indexList = list(range(len(genome)))
    sortedIndex = merge_sort(indexList)
    return sortedIndex
    

def bwt(t):
    """ Given T, returns BWT(T), by way of the suffix array. """
    return ''.join([t[si-1] for si in suffixArray(t)])

# 
bwtgenome = bwt(genome)

# bwtgenome = 'smnpbnnaaaaa$a'
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
# freq = {letter:0 for letter in alphabets}
count = {letter:[0]*(genomeLen+1) for letter in alphabets}
for i in range(genomeLen):
    for letter in alphabets:
        count[letter][i+1] = count[letter][i]
    count[bwtgenome[i]][i+1] += 1
    
#print(count)

def bwtMatch(pattern):
    top = 0
    bottom = genomeLen - 1
    initialPattern = pattern

    while top <= bottom:
        if len(pattern) > 0:
            symbol = pattern[-1]
            pattern = pattern[0:-1]
            if count[symbol][bottom+1] - count[symbol][top] > 0:
                top = firstOccurence[symbol] + count[symbol][top]
                bottom = firstOccurence[symbol] + count[symbol][bottom+1] - 1
            else:
                return []
        else:
            #if top >= 5792 and bottom <= 5792:
            #    print(initialPattern)
            return [x for x in range(top, bottom+1)]


matchBwtIndexes = [bwtMatch(p) for p in patterns]
#print(matchBwtIndexes)
# flatten the list
matchBwtIndexes = list(itertools.chain.from_iterable(matchBwtIndexes))
# take only unique matches
matchBwtIndexes = set(matchBwtIndexes)
#print(matchBwtIndexes)

# convert to original index
sarr = suffixArray(genome)
matchBwtIndexes = [sarr[idx] for idx in matchBwtIndexes]
matchBwtIndexes = sorted(matchBwtIndexes)


## output
with open(outputFile, "w") as f:
    f.writelines(' '.join([str(n) for n in matchBwtIndexes]) + '\n')

print('done.')
