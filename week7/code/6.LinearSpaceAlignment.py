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
inputFile = './data/6.LinearSpaceAlignment-2.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_250_14.txt'
outputFile = './results/6.LinearSpaceAlignment.txt'
scoreFile = './data/BLOSUM62.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

protein1 = inputs[0].strip()
protein2 = inputs[1].strip()

#protein2 = "MEASNLY"

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


def midEdge(top, bottom, left, right):
    global protein1, protein2
    global prevScores, curScores
    global blosum, colIndex
    
    middle = math.floor((left+right)/2)

    ## initialize score for source to middle
    prevScores[top] = 0
    for i in range(top+1,bottom+1):
        prevScores[i] = prevScores[i-1] - sigma

    ## from source to middle
    for j in range(left+1, middle+1):
        curScores[top] = prevScores[top] - sigma
        for i in range(top+1,bottom+1):
            topScore = curScores[i-1] - sigma
            leftScore = prevScores[i] - sigma
            diagScore = prevScores[i-1] + blosum[colIndex[protein1[i-1]]][colIndex[protein2[j-1]]]
            curScores[i] = max(topScore, leftScore, diagScore)

        # swap prev and cur scores
        prevScores, curScores = curScores, prevScores

    fromSourceScores = prevScores.copy()


    ## initialize score for sink to middle
    prevScores[bottom] = 0
    for i in range(bottom-1,top-1,-1):
        prevScores[i] = prevScores[i+1] - sigma

    ## from sink to middle
    for j in range(right-1,middle-1, -1):
        curScores[bottom] = prevScores[bottom] - sigma
        for i in range(bottom-1,top-1,-1):
            bottomScore = curScores[i+1] - sigma
            rightScore = prevScores[i] - sigma
            diagScore = prevScores[i+1] + blosum[colIndex[protein1[i]]][colIndex[protein2[j]]]
            curScores[i] = max(bottomScore, rightScore, diagScore)

        # swap prev and cur scores
        prevScores, curScores = curScores, prevScores

    middleToSinkScores = prevScores.copy()
    nextToSinkScore = curScores.copy()


    ## find middle node
    scores = [fromSourceScores[i]+middleToSinkScores[i] for i in range(top, bottom+1)]
    maxScore = max(scores)
    middleNode = [top+scores.index(maxScore), middle]
    #print(middleNode)

    ## find middle edge
    nextNode = [0,0]
    midRow = middleNode[0]
    midCol = middleNode[1]
    if(midRow==bottom):
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
                
    #print(middleNode, nextNode)

    return [middleNode, nextNode, maxScore]

alignment1 = []
alignment2 = []
def linearSpaceAlignment(top, bottom, left, right):
    global protein1, protein2
    global prevScores, curScores
    global blosum, colIndex
    global alignment1, alignment2
    
    if left == right:
        a1 = [protein1[i] for i in range(top,bottom)]
        a2 = ['-']*(bottom-top)
        alignment1.extend(a1)
        alignment2.extend(a2)
        return -sigma*(bottom-top)

    if top == bottom:
        a1 = ['-']*(right-left)
        a2 = [protein2[i] for i in range(left,right)]
        alignment1.extend(a1)
        alignment2.extend(a2)
        return -sigma*(right-left)

    me = midEdge(top, bottom, left, right)
    middle = math.floor((left+right)/2)
    linearSpaceAlignment(top, me[0][0], left, middle)

    if me[0][0]==me[1][0] and me[0][1]+1==me[1][1]:     #right
        alignment1.append('-')
        alignment2.append(protein2[me[0][1]])
    elif me[0][0]+1==me[1][0] and me[0][1]==me[1][1]:     #bottom
        alignment1.append(protein1[me[0][0]])
        alignment2.append('-')
    else:
        alignment1.append(protein1[me[0][0]])
        alignment2.append(protein2[me[0][1]])

    linearSpaceAlignment(me[1][0], bottom, me[1][1], right)
    return me[2]


alignmentScore = linearSpaceAlignment(0, len1, 0, len2)
print(alignmentScore)

#me = midEdge(0, len1, 0, len2)
#print(me)


    



#alignmentScore = 0
#for i in range(len(alignment1)):
#    if alignment1[i]=='-' or alignment2[i]=='-' :
#        alignmentScore = alignmentScore - sigma
#        #print(str(-sigma), end=" ")
#    else:
#        alignmentScore = alignmentScore + blosum[colIndex[alignment1[i]]][colIndex[alignment2[i]]]
#        #print(str(blosum[colIndex[alignment1[i]]][colIndex[alignment2[i]]]), end=" ")
#
#print(alignmentScore)
#print(alignment1)
#print(alignment2)


# output
with open(outputFile, "w") as f:
    f.writelines(str(alignmentScore) + "\n")
    f.writelines("".join(alignment1) + "\n")
    f.writelines("".join(alignment2) + "\n")
    
