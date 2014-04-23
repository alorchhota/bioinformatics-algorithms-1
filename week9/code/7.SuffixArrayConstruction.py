import os
import csv
import sys
import re
import importlib
import networkx as nx

# settings
#curDir = 'E:/GitHub/bioinformatics-algorithms-1/week9'
curDir = 'D:/GitHub/bioinformatics-algorithms-1/week9'
inputFile = './data/7.SuffixArrayConstruction-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_310_2.txt'
outputFile = './results/7.SuffixArrayConstruction.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome = inputs[0].strip()
#genome = 'GFEDCBA'
#genome = genome[0:50000]
genomeLen = len(genome)

suffix  = lambda idx: genome[idx:]
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

#tree = suffixTree(genome)

#indexList = list(range(genomeLen))
#indexList.sort(key=lambda x: genome[x:])

#genome = 'AACGATAGCGGTAGA$'

    

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

indexList = list(range(genomeLen))

sortedList = merge_sort(indexList)
#print(sortedList)
## output
with open(outputFile, "w") as f:
    f.writelines(', '.join([str(i) for i in sortedList]))

print('done.')
