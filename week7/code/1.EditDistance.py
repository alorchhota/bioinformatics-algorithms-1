import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
inputFile = './data/1.EditDistance-1.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_248_3.txt'
outputFile = './results/1.EditDistance.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

protein1 = inputs[0].strip()
protein2 = inputs[1].strip()

sigma = 1


## create scores and backtrack arrays
len1 = len(protein1)
len2 = len(protein2)
scores = [[0]*(len2+1) for x in range(len1+1)]

## put first row and first column of scores and backtrack arrays
for i in range(1,len1+1):
    scores[i][0] = scores[i-1][0] + sigma
    
for j in range(1,len2+1):
    scores[0][j] = scores[0][j-1] + sigma

## update scores in greedy approach
for i in range(1,len1+1):
    for j in range(1,len2+1):
        topScore = scores[i-1][j] + sigma
        leftScore = scores[i][j-1] + sigma
        diagScore = scores[i-1][j-1] + int(protein1[i-1]!=protein2[j-1])
        candidateScores = [diagScore, topScore, leftScore]
        scores[i][j] = min(candidateScores)

## max score
minScore = scores[len1][len2]

print(minScore)

# output
with open(outputFile, "w") as f:
    f.writelines(str(minScore) + "\n")
