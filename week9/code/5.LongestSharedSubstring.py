import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
curDir = 'E:/GitHub/bioinformatics-algorithms-1/week9'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week9'
inputFile = './data/5.LongestSharedSubstring-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_296_5.txt'
outputFile = './results/5.LongestSharedSubstring.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome1 = inputs[0].strip() + '$'
genome2 = inputs[1].strip() + '$'

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

## build two suffix tree for two genomes
g1 = suffixTree(genome1)
g2 = suffixTree(genome2)

## function to find edges with prefix from a node
neighborsWithLabelPrefix = lambda genome, g, node, prefix:\
                           [e[1] for e in g.edges_iter(node, data=True) \
                            if genome[e[2]['labelIdx'][0]] == prefix]

## function to find edge label from suffix tree
edgeLabelInSuffixTree = lambda genome, g, startNode, endNode:\
                        genome[g.edge[startNode][endNode]['labelIdx'][0]:\
                               g.edge[startNode][endNode]['labelIdx'][1]]

def longestSharedSubstring(root1, root2):
    longestSubstr = ''

    ## create matched edge pairs
    n1 = nx.neighbors(g1, root1)
    if len(n1) == 0:
        return longestSubstr

    edges1 = [(edgeLabelInSuffixTree(genome1, g1, root1, node), node) for node in n1]
    edgePairs = []
    for e in edges1:
        n2 = neighborsWithLabelPrefix(genome2, g2, root2, e[0][0])
        if len(n2)>0:
            e2 = (edgeLabelInSuffixTree(genome2, g2, root2, n2[0]), n2[0])
            edgePairs += [(e,e2)]

    if len(edgePairs) == 0:
        return longestSubstr

    ## traverse each edge pairs and update longest substr
    for ep in edgePairs:

        ## find substr in each pair
        substr = ''
        cur1 = root1
        cur2 = root2
        next1 = ep[0][1]
        next2 = ep[1][1]
        edge1 = ep[0][0]
        edge2 = ep[1][0]

        while True:
            if edge1 == edge2:
                substr = edge1 + longestSharedSubstring(next1, next2)
                break

            # update substr with prefix
            prefix = longestCommonPrefix(edge1, edge2)
            substr += prefix

            if len(edge1) < len(edge2):
                edge2 = edge2[len(prefix):]
                cur1 = next1
                next1 = neighborsWithLabelPrefix(genome1, g1, cur1, edge2[0])
                if len(next1) == 0:
                    break
                next1 = next1[0]
                edge1 = edgeLabelInSuffixTree(genome1, g1, cur1, next1)
            else:
                edge1 = edge1[len(prefix):]
                cur2 = next2
                next2 = neighborsWithLabelPrefix(genome2, g2, cur2, edge1[0])
                if len(next2) == 0:
                    break
                next2 = next2[0]
                edge2 = edgeLabelInSuffixTree(genome2, g2, cur2, next2)

        # update longest substring
        if len(substr) > len(longestSubstr):
            longestSubstr = substr

    return longestSubstr
            
                
                
    

## find longest shared sustring
lss = longestSharedSubstring(1,1)
print(lss)

## output
with open(outputFile, "w") as f:
    f.writelines(lss)

print('done.')
