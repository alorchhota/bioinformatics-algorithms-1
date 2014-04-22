import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/GitHub/bioinformatics-algorithms-1/week9'
curDir = 'D:/GitHub/bioinformatics-algorithms-1/week9'
inputFile = './data/6.ShortestNonSharedSubstring-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_296_6 (1).txt'
outputFile = './results/6.ShortestNonSharedSubstring.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome1 = inputs[0].strip() 
genome2 = inputs[1].strip() + '$'

genomeLen1 = len(genome1)
genomeLen2 = len(genome2)

## function to find longest common prefix
def longestCommonPrefix(str1, str2):
    n = min([len(str1), len(str2)])
    i = 0
    while i < n and str1[i]==str2[i]:
        i += 1

    prefix = str1[0:i]
    return prefix

## function to build a suffix tree from a genome
def suffixTree(genome):
    ## build suffix tree
    g = nx.DiGraph()
    g.add_node(1)     # add root with id 1

    # two required function
    neighborsWithLabelPrefix = lambda node, prefix: [e[1] for e in g.edges_iter(node, data=True) if genome[e[2]['labelIdx'][0]] == prefix]
    getNewNode = lambda : len(g.nodes())+1

    #print(longestCommonPrefix('abc','ab'))

    genomeLen = len(genome)
    for idx in range(genomeLen):
        # traverse as long as pattern matches
        curNode = 1
        i = idx
        while(i < genomeLen):
            # find the edge with the first prefix character
            nextNode = neighborsWithLabelPrefix(curNode, genome[i])

            # if there is no edge with the first prefix character,
            # it must be a new edge with the rest of the string.
            if len(nextNode) == 0:
                newNode = getNewNode()
                g.add_edge(curNode, newNode, {'labelIdx':[i,genomeLen]})
                g.node[newNode]['startIdx'] = idx
                break

            # get the edge label
            nextNode = nextNode[0]
            edgeLabelIndices = g.edge[curNode][nextNode]['labelIdx']
            edgeLabel = genome[edgeLabelIndices[0]:edgeLabelIndices[1]]
            edgeLabelLen = len(edgeLabel) 

            # if the rest of the string starts with edgeLabel,
            # move to the next node
            if genome[i:i+edgeLabelLen] == edgeLabel:
                curNode = nextNode
                i += edgeLabelLen
            else:
                # edgeLabel matches partially
                prefix = longestCommonPrefix(genome[i:i+edgeLabelLen], edgeLabel)
                prefixLen = len(prefix)

                # create two new node, one intermediate, another for unmatched string
                intermediateNode = getNewNode()
                unmatchedNode = intermediateNode + 1

                # remove existing edge from curNode to nextNode
                g.remove_edge(curNode, nextNode)
                
                # add edge from curNode to intermediateNode
                g.add_edge(curNode, intermediateNode, {'labelIdx':(edgeLabelIndices[0],edgeLabelIndices[0]+prefixLen)})

                # add edge from intermediateNode to nextNode
                g.add_edge(intermediateNode, nextNode, {'labelIdx':(edgeLabelIndices[0]+prefixLen, edgeLabelIndices[1])})

                # add edge from intermediateNode to unmatchedNode
                g.add_edge(intermediateNode, unmatchedNode, {'labelIdx':(i+prefixLen, genomeLen)})
                g.node[unmatchedNode]['startIdx'] = idx
                
                break

    return g

## build suffix tree of the first genome
g2 = suffixTree(genome2)

## function to find edges with prefix from a node
neighborsWithLabelPrefix = lambda genome, g, node, prefix:\
                           [e[1] for e in g.edges_iter(node, data=True) \
                            if genome[e[2]['labelIdx'][0]] == prefix]

## function to find edge label from suffix tree
edgeLabelInSuffixTree = lambda genome, g, startNode, endNode:\
                        genome[g.edge[startNode][endNode]['labelIdx'][0]:\
                               g.edge[startNode][endNode]['labelIdx'][1]]


## define function to find shorted non-shared node from a position
def shortestNonSharedSubstrLengthFrom(idx, maxLength):
    cur2 = 1    # start at the root of g2
    i = idx
    lim = min([idx+maxLength, genomeLen1])
    #print("idx="+str(idx))
    while i < lim:
        #print(cur2, genome1[i])
        nextNodes = neighborsWithLabelPrefix(genome2, g2, cur2, genome1[i])
        if len(nextNodes) == 0:
            # did not match
            return i-idx+1

        nextNode = nextNodes[0]
        #print(nextNode)
        edge = edgeLabelInSuffixTree(genome2, g2, cur2, nextNode)
        edgeLen = len(edge)
        if edge == genome1[i:i+edgeLen]:
            cur2 = nextNode
            i += edgeLen
        else:
            prefix = longestCommonPrefix(edge, genome1[i:i+edgeLen])
            prefixLen = len(prefix)
            i += prefixLen

            if i >= lim:
                # matched fully
                return float('inf')

            return i+1-idx

    # matched fully
    return float('inf')

## find shortest non shared node
shortestSubstrStartIndex = -1
shortestLength = genomeLen1

#edges = [(e[0],e[1], edgeLabelInSuffixTree(genome2,g2,e[0],e[1])) for e in g2.edges_iter(data=True)]
#print(edges)

for i in range(genomeLen1):
    length = shortestNonSharedSubstrLengthFrom(i, genomeLen1)
    #print(length)
    if i % 500 == 0:
        print(i, shortestLength)
    if length < shortestLength:
        shortestLength = length
        shortestSubstrStartIndex = i
    
shortestNonSharedSubstr = genome1[shortestSubstrStartIndex:shortestSubstrStartIndex+shortestLength]
print(shortestNonSharedSubstr)

## output
with open(outputFile, "w") as f:
    f.writelines(shortestNonSharedSubstr)

print('done.')
