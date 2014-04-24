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
inputFile = './data/5.PartialSuffixArray-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_302_3.txt'
outputFile = './results/5.PartialSuffixArray.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome = inputs[0].strip()
k = int(inputs[1].strip())

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
    
def partialSuffixArray(genome, k):
    sarr = suffixArray(genome)
    psarr = [(i, sarr[i]) for i in range(len(sarr)) if sarr[i]%k==0]
    return psarr


psarr = partialSuffixArray(genome, k)
print('\n'.join([','.join([str(i) for i in item])  for item in psarr]))

## output
with open(outputFile, "w") as f:
    f.writelines('\n'.join([','.join([str(i) for i in item])  for item in psarr]))

print('done.')
