import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
inputFile = './data/4.affine_gap-3.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_249_8.txt'
outputFile = './results/4.affine_gap.txt'
scoreFile = './data/BLOSUM62.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

protein1 = inputs[0].strip()
protein2 = inputs[1].strip()

sigma = 11
epsilon = 1

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
# [0][][]:lower, [1][][]:mid, [2][][]:upper
scores = [[[0]*(len2+1) for x in range(len1+1)],\
         [[0]*(len2+1) for x in range(len1+1)],\
         [[0]*(len2+1) for x in range(len1+1)]]

# 0: top(lower)     | diag(mid)         | left(upper)
# 1: from-mid(lower)| from-lower(mid)   | from-mid(upper)
# 2:  -             | from-upper(mid)   |
backtrack = [[[0]*(len2+1) for x in range(len1+1)],\
             [[0]*(len2+1) for x in range(len1+1)],\
             [[0]*(len2+1) for x in range(len1+1)]]


## put first row and first column of scores and backtrack arrays
## first row and cols are negative infinity,
## except the first entry of the mid scores
for i in range(0,len1+1):
    scores[0][i][0] = float('-inf')
    scores[1][i][0] = float('-inf')
    scores[2][i][0] = float('-inf')
for j in range(0,len2+1):
    scores[0][0][j] = float('-inf')
    scores[1][0][j] = float('-inf')
    scores[2][0][j] = float('-inf')

scores[1][0][0] = 0

## update scores in greedy approach
for i in range(1,len1+1):
    for j in range(1,len2+1):
        # lower score
        topScore = scores[0][i-1][j] - epsilon
        fromMidScore = scores[1][i-1][j] - sigma
        candidateScores = [topScore, fromMidScore]
        scores[0][i][j] = max(candidateScores)
        backtrack[0][i][j] = candidateScores.index(scores[0][i][j])
        
        # upper score
        leftScore = scores[2][i][j-1] - epsilon
        fromMidScore = scores[1][i][j-1] - sigma
        candidateScores = [leftScore, fromMidScore]
        scores[2][i][j] = max(candidateScores)
        backtrack[2][i][j] = candidateScores.index(scores[2][i][j])

        # mid score
        diagScore = scores[1][i-1][j-1] + blosum[colIndex[protein1[i-1]]][colIndex[protein2[j-1]]]
        fromLowerScore = scores[0][i][j]
        fromUpperScore = scores[2][i][j]
        candidateScores = [diagScore, fromLowerScore, fromUpperScore]
        scores[1][i][j] = max(candidateScores)
        backtrack[1][i][j] = candidateScores.index(scores[1][i][j])
        

## max score is found in the mid level
maxScore = scores[1][len1][len2]

## backtrack and find alignment
alignment1 = []
alignment2 = []

# 0: top(lower)     | diag(mid)         | left(upper)
# 1: from-mid(lower)| from-lower(mid)   | from-mid(upper)
# 2:  -             | from-upper(mid)   |

k = 1       # which level of scoring lower/mid/high
i = len1    # row index
j = len2    # col index
while (k==1 and i==0 and j==0) != True:
    if k == 0:      # in lower scores
        alignment1.append(protein1[i-1])
        alignment2.append('-')
        if backtrack[k][i][j] == 0:     # top
            i = i-1
        else:                           # from-mid
            i = i-1
            k = 1
    elif k == 1:    # in mid scores
        if backtrack[k][i][j] == 0:     # diag
            alignment1.append(protein1[i-1])
            alignment2.append(protein2[j-1])
            i = i-1
            j = j-1
        elif backtrack[k][i][j] == 1:   # from-lower
            k = 0
        else:                           # from-upper
            k = 2
    else:           # in upper scores
        alignment1.append('-')
        alignment2.append(protein2[j-1])
        if backtrack[k][i][j] == 0:     # left
            j = j-1
        else:                           # from-mid
            j = j-1
            k = 1

    
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

