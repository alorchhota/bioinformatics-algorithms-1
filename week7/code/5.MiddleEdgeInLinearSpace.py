import os
import csv
import sys
import re
import importlib
#import networkx as nx
import math

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
inputFile = './data/5.MiddleEdgeInLinearSpace-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_250_12.txt'
outputFile = './results/5.MiddleEdgeInLinearSpace.txt'
scoreFile = './data/BLOSUM62.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

protein1 = inputs[0].strip()
protein2 = inputs[1].strip()

sigma = 5

## read blosum score matrix
## colIndex is map from amino acid to its index
## blosum is a 2D matrix
with open(scoreFile) as f:
    lines = f.readlines()
    cols = re.split("\\s*", lines[0].strip())
    colIndex = {cols[i]:i for i in range(0,len(cols))}
    blosum = [[int(i) for i in re.split("\\s*", row.strip())[1:]] for row in lines[1:]]


## create scores and backtrack arrays
len1 = len(protein1)
len2 = len(protein2)

## scores storage and initialization
prevScores = [0]*(len1+1)   # scores in previous col
curScores = [0]*(len1+1)    # # scores in current col

prevScores[0] = 0
for i in range(1,len1+1):
    prevScores[i] = prevScores[i-1] - sigma

## middle edge
middle = math.floor(len2/2)

## from source to middle
for j in range(1,middle+1):
    curScores[0] = prevScores[0] - sigma
    for i in range(1,len1+1):
        topScore = curScores[i-1] - sigma
        leftScore = prevScores[i] - sigma
        diagScore = prevScores[i-1] + blosum[colIndex[protein1[i-1]]][colIndex[protein2[j-1]]]
        curScores[i] = max(topScore, leftScore, diagScore)

    # swap prev and cur scores
    prevScores, curScores = curScores, prevScores

fromSourceScores = prevScores.copy()


## from sink to middle
prevScores[len1] = 0
for i in range(len1-1,-1,-1):
    prevScores[i] = prevScores[i+1] - sigma
    
for j in range(len2-1,middle-1, -1):
    curScores[len1] = prevScores[len1] - sigma
    for i in range(len1-1,-1,-1):
        bottomScore = curScores[i+1] - sigma
        rightScore = prevScores[i] - sigma
        diagScore = prevScores[i+1] + blosum[colIndex[protein1[i]]][colIndex[protein2[j]]]
        curScores[i] = max(bottomScore, rightScore, diagScore)

    # swap prev and cur scores
    prevScores, curScores = curScores, prevScores

middleToSinkScores = prevScores.copy()
nextToSinkScore = curScores.copy()


## find middle node
scores = [fromSourceScores[i]+middleToSinkScores[i] for i in range(len1+1)]
middleNode = [scores.index(max(scores)), middle]
#print(max(scores))

## find middle edge
nextNode = [0,0]
midRow = middleNode[0]
midCol = middleNode[1]
if(midRow==len1):
    nextNode = [midRow, midCol+1] 
else:
    rightScore = - sigma + nextToSinkScore[midRow]
    bottomScore = - sigma + middleToSinkScores[midRow+1]
    diagScore = blosum[colIndex[protein1[midRow]]][colIndex[protein2[midCol]]] + nextToSinkScore[midRow+1]
    nextScore = max(rightScore, bottomScore, diagScore)
    if nextScore == rightScore:
        nextNode = [midRow, midCol+1]
    elif nextScore == bottomScore:
        nextNode = [midRow+1, midCol]
    else:
        nextNode = [midRow+1, midCol+1]
            
print(middleNode, nextNode)

# output
with open(outputFile, "w") as f:
    f.writelines( '(' + str(middleNode[0]) + ',' + str(middleNode[1]) + ') (' + str(nextNode[0]) + ',' + str(nextNode[1]) + ')')

