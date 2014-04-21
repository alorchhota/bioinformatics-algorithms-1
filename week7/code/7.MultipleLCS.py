import os
import csv
import sys
import re
import importlib
#import networkx as nx

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week7'
inputFile = './data/7.MultipleLCS-1.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_251_5.txt'
outputFile = './results/7.MultipleLCS.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

protein1 = inputs[0].strip()
protein2 = inputs[1].strip()
protein3 = inputs[2].strip()


## create scores and backtrack arrays
len1 = len(protein1)
len2 = len(protein2)
len3 = len(protein3)
scores = [[[0]*(len3+1) for x in range(len2+1)] for y in range(len1+1)]
backtrack = [[[0]*(len3+1) for x in range(len2+1)] for y in range(len1+1)]


## put first row and first column of scores and backtrack arrays
scores[0][0][0] = 0
for i in range(1,len1+1):
    scores[i][0][0] = scores[i-1][0][0] - 1
    backtrack[i][0][0] = 0
    
for j in range(1,len2+1):
    scores[0][j][0] = scores[0][j-1][0] - 1
    backtrack[0][j][0] = 1

for k in range(1,len3+1):
    scores[0][0][k] = scores[0][0][k-1] - 1
    backtrack[0][0][k] = 2
    


## update scores in greedy approach
for i in range(1,len1+1):
    for j in range(1,len2+1):
        for k in range(1, len3+1):
            s0 = scores[i-1][j][k]
            s1 = scores[i][j-1][k]
            s2 = scores[i][j][k-1]
            s3 = scores[i-1][j-1][k]
            s4 = scores[i-1][j][k-1]
            s5 = scores[i][j-1][k-1]
            s6 = scores[i-1][j-1][k-1] + int(protein1[i-1]==protein2[j-1]==protein3[k-1])
            candidateScores = [s0, s1, s2, s3, s4, s5, s6]
            scores[i][j][k] = max(candidateScores)
            backtrack[i][j][k] = candidateScores.index(scores[i][j][k])

## max score node is sink score
maxScore = scores[len1][len2][len3]
print(str(maxScore))


## backtrack and find alignment
alignment1 = []
alignment2 = []
alignment3 = []

i = len1
j = len2
k = len3
while i!=0 or j!=0 or k!=0:
    if backtrack[i][j][k] == 0:
        alignment1.append(protein1[i-1])
        alignment2.append('-')
        alignment3.append('-')
        i = i-1
    elif backtrack[i][j][k] == 1:
        alignment1.append('-')
        alignment2.append(protein2[j-1])
        alignment3.append('-')
        j = j-1
    elif backtrack[i][j][k] == 2:
        alignment1.append('-')
        alignment2.append('-')
        alignment3.append(protein3[k-1])
        k = k-1
    elif backtrack[i][j][k] == 3:
        alignment1.append(protein1[i-1])
        alignment2.append(protein2[j-1])
        alignment3.append('-')
        i = i-1
        j = j-1
    elif backtrack[i][j][k] == 4:
        alignment1.append(protein1[i-1])
        alignment2.append('-')
        alignment3.append(protein3[k-1])
        i = i-1
        k = k-1
    elif backtrack[i][j][k] == 5:
        alignment1.append('-')
        alignment2.append(protein2[j-1])
        alignment3.append(protein3[k-1])
        j = j-1
        k = k-1
    else:
        alignment1.append(protein1[i-1])
        alignment2.append(protein2[j-1])
        alignment3.append(protein3[k-1])
        i = i-1
        j = j-1
        k = k-1 
    
alignment1.reverse()
alignment2.reverse()
alignment3.reverse()

print(maxScore)
print(alignment1)
print(alignment2)
print(alignment3)

# output
with open(outputFile, "w") as f:
    f.writelines(str(maxScore) + "\n")
    f.writelines("".join(alignment1) + "\n")
    f.writelines("".join(alignment2) + "\n")
    f.writelines("".join(alignment3) + "\n")

