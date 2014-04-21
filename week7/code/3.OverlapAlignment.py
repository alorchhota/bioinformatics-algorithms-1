import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
inputFile = './data/3.overlap_alignment-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_248_7.txt'
outputFile = './results/3.overlap_alignment.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

protein1 = inputs[0].strip()
protein2 = inputs[1].strip()

sigma = 2



## create scores and backtrack arrays
len1 = len(protein1)
len2 = len(protein2)
scores = [[0]*(len2+1) for x in range(len1+1)]
backtrack = [[0]*(len2+1) for x in range(len1+1)]     # 0: top, 1:left, 2:diag, 3:start


## put first row and first column of scores and backtrack arrays
for i in range(1,len1+1):
    scores[i][0] = 0
    backtrack[i][0] = 3     # start
    
for j in range(1,len2+1):
    scores[0][j] = scores[0][j-1] - sigma
    backtrack[0][j] = 1


## update scores in greedy approach
for i in range(1,len1+1):
    for j in range(1,len2+1):
        topScore = scores[i-1][j] - sigma
        leftScore = scores[i][j-1] - sigma
        if protein1[i-1] == protein2[j-1]:
            diagExtraScore = 1
        else:
            diagExtraScore = -2
        diagScore = scores[i-1][j-1] + diagExtraScore
        candidateScores = [topScore, leftScore, diagScore]
        scores[i][j] = max(candidateScores)
        backtrack[i][j] = candidateScores.index(scores[i][j])

## max score node is the max in the last row where 1st protein ends
maxScore = float("-inf")
maxj = -1
for j in range(len2+1):
    if scores[len1][j] > maxScore:
        maxScore = scores[len1][j]
        maxj = j


## backtrack and find alignment
alignment1 = []
alignment2 = []

i = len1
j = maxj
while i!=0 or j!=0:
    if backtrack[i][j] == 0:        # top
        alignment1.append(protein1[i-1])
        alignment2.append('-')
        i = i-1
    elif backtrack[i][j] == 1:      # left
        alignment1.append('-')
        alignment2.append(protein2[j-1])
        j = j-1
    elif backtrack[i][j] == 2:      # diag
        alignment1.append(protein1[i-1])
        alignment2.append(protein2[j-1])
        i = i-1
        j = j-1
    else:                           # start
        i = 0
        j = 0
    
alignment1.reverse()
alignment2.reverse()

print(maxScore)
#print(alignment1)
#print(alignment2)

# output
with open(outputFile, "w") as f:
    f.writelines(str(maxScore) + "\n")
    f.writelines("".join(alignment1) + "\n")
    f.writelines("".join(alignment2) + "\n")

