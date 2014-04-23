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
inputFile = './data/1.bwt-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_297_4.txt'
outputFile = './results/1.bwt.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome = inputs[0].strip()
genomeLen = len(genome)


def suffixArray_v1(s):
    satups = sorted([(s[i:], i) for i in range(len(s))])
    # Extract and return just the offsets
    return map(lambda x: x[1], satups)

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

bwtGenome = bwt(genome)
#print(bwtGenome)

## output
with open(outputFile, "w") as f:
    f.writelines(bwtGenome + '\n')

print('done.')
