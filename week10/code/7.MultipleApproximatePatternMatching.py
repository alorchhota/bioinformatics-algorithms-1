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
inputFile = './data/7.MultipleApproximatePatternMatching-2.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_304_6.txt'
outputFile = './results/7.MultipleApproximatePatternMatching.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome = inputs[0].strip() + '$'
patterns = [p for p in re.split('\\s+', inputs[1].strip())]
d = int(inputs[2].strip())

#print(genome)
#print (patterns)
#print(d)

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

def hammingDistance(s1, s2):
    if len(s1) != len(s2):
        return float("inf")
    return len([i for i in range(len(s1)) if s1[i]!=s2[i]])

    

def bwtPartialMatch(pattern, d, sarr):
    # divide the pattern into d+1 parts
    patternLen = len(pattern)
    divLen = int(patternLen/(d+1))
    extra = patternLen % (d+1)
    divLen = [divLen]*(d+1-extra) + [divLen+1]*extra
    divIndex = list(itertools.accumulate([0]+divLen))

    # find a pattern that has exact match
    smallPatterns = [pattern[divIndex[pi]:divIndex[pi+1]] for pi in range(d+1)]
    smallMatchIndexes = [bwtMatch(p) for p in smallPatterns]

    # extend the non empty patterns
    nonEmptyPatternMatchIndexes = [pi for pi in range(d+1) if len(smallMatchIndexes[pi]) > 0]
    allPartialMatchIndexes = set()
    for pi in range(d+1):
        if len(smallMatchIndexes[pi]) > 0:
            # if there are exact matches
            exactMatchIndexes = smallMatchIndexes[pi]
            extendedStartIndexes = [sarr[mi] - divIndex[pi] \
                                    for mi in exactMatchIndexes\
                                    if sarr[mi] - divIndex[pi] >= 0]
            
            partialMatchIndexes = [si for si in extendedStartIndexes \
                                   if hammingDistance(genome[si:si+patternLen], pattern) <= d]

            for i in partialMatchIndexes:
                allPartialMatchIndexes.add(i)
            
        
    return allPartialMatchIndexes

            
                    
    
    
    
sarr = suffixArray(genome)    
matchBwtIndexes = [bwtPartialMatch(p, d, sarr) for p in patterns]
#print(matchBwtIndexes)
# flatten the list
matchBwtIndexes = list(itertools.chain.from_iterable(matchBwtIndexes))
# take only unique matches
#matchBwtIndexes = set(matchBwtIndexes)
matchBwtIndexes = sorted(matchBwtIndexes)
#print(matchBwtIndexes)

# convert to original index

#matchBwtIndexes = [sarr[idx] for idx in matchBwtIndexes]
#matchBwtIndexes = sorted(matchBwtIndexes)


## output
with open(outputFile, "w") as f:
    f.writelines(' '.join([str(n) for n in matchBwtIndexes]) + '\n')

print('done.')
