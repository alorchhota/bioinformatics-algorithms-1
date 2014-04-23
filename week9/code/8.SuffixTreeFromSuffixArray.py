import os
import csv
import sys
import re
import importlib
import networkx as nx
import itertools

# settings
curDir = 'E:/GitHub/bioinformatics-algorithms-1/week9'
#curDir = 'D:/GitHub/bioinformatics-algorithms-1/week9'
inputFile = './data/8.SuffixTreeFromSuffixArray-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_310_5.txt'
outputFile = './results/8.SuffixTreeFromSuffixArray.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

genome1 = inputs[0].strip()
suffixArray1 = [int(i) for i in re.split(',\\s*', inputs[1].strip())]
lcp1 = [int(i) for i in re.split(',\\s*', inputs[2].strip())]
# print(genome, suffixArray, lcp)

#genomeLen = len(genome)

## settings
edgeLabelAttr = 'labelinfo'
nodeMatchAttr = 'matchindex'    # keeps the start index where match started


###### algorithm ######
## 0.keep a storages for
##     1.unconnected-nodes
##     2.label from previous unconnected nodes to the current one.
## 1.create a root node
## 2.move through suffix array lexicographically.
##     1.if prevLCP < nextLCP (i.e. some portion matches in forward direction)
##        1.create a new node
##        2.add curNode to unconnected nodes with forward matched label.
##        3.move to newly connected node (curNode = new node)
##        4.if there is unmatched tail,\
##          add an edge from the cur node with that label.
##    2. if prevLCP == nextLCP (no forward match, no reverse match)
##        1.if there is unmatched tail,\
##          add an edge from the cur node with that label.
##    3. if prevLCP > nextLCP (some reverse match)
##        1. if there is unmatched tail,\
##          add an edge from the cur node with that label.
##        2. You have to move backward in unconnected nodes and labels\
##           until reverse match is consumed by saved label.\
##           while |lastLabel| <= |unconsumedMatch|:
##               add an edge from lastNode to curNode with lastLabel
##               move to lasNode (curNode = lastNode)
##               reduce unconsumedMatch by |lastLabel|
##           if |unconsumedMatch| > 0: # this is an intermediate node
##               create an intermediateNode
##               get the lastNode with labels
##               create an edge from intermediateNode to curNode with unconsumedMatch
##               move to intermediateNode
##               change the label of lastNode to lastLabel - unconsumedMatch
##               

def suffixArrayToSuffixTree (genome, suffixArray, lcp):

    g = nx.DiGraph()

    ## add an extra zero to the end of lcp for programming purpose
    lcp.append(0)

    genomeLen = len(genome)
    unconnectedNodes = []
    unconnectedLabels = []      # (startIndexInGenome, length)

    ## utility functions
    getNewNode = lambda : len(g.nodes())+1
    lengthOfSuffix = lambda idx: genomeLen - idx


    ## add root
    g.add_node(1)
    #unconnectedNodes.append(1)
    #unconnectedLabels.append((0,0))

    ## move through suffix array lexicographically.
    curNode = 1 # root
    for idx in range(len(suffixArray)):
        sidx = suffixArray[idx]
        if lcp[idx] < lcp[idx+1]:   # some portion matches in forward direction
            # create a node, add to unconnected nodes with forward matched label.
            newNode = getNewNode()
            g.add_node(newNode)
            unconnectedNodes.append(curNode)
            matchedLabel = (sidx+lcp[idx], lcp[idx+1]-lcp[idx])
            unconnectedLabels.append(matchedLabel)

            # move to newly connected node
            curNode = newNode
            
        # if there is unmatched tail,\
        # add an edge from the cur node with that label.
        unmatchedStartIdx = max([lcp[idx], lcp[idx+1]])
        suffixLen = lengthOfSuffix(sidx)
        if unmatchedStartIdx < suffixLen:
            newNode = getNewNode()
            g.add_node(newNode, {nodeMatchAttr:sidx})
            g.add_edge(curNode, newNode, {edgeLabelAttr:(sidx+unmatchedStartIdx, suffixLen-unmatchedStartIdx)})
            

        if lcp[idx] > lcp[idx+1]:    # some reverse match
            # move backward in unconnected nodes and labels
            # until reverse match is consumed by saved label.
            reverseMatchLen = lcp[idx]-lcp[idx+1]
            while len(unconnectedLabels)>0 and unconnectedLabels[-1][1] <= reverseMatchLen:
                #add an edge from lastNode to curNode with lastLabel
                lastNode = unconnectedNodes.pop()
                lastLabel = unconnectedLabels.pop()
                g.add_edge(lastNode, curNode, {edgeLabelAttr:lastLabel})
                
                # move to lasNode (curNode = lastNode)
                curNode = lastNode
                
                # reduce unconsumedMatch by |lastLabel|
                reverseMatchLen -= lastLabel[1]


            if reverseMatchLen > 0:
                # create an intermediateNode
                interNode = getNewNode()
                g.add_node(interNode)
                
                #get the lastNode with labels
                lastNode = unconnectedNodes[-1]
                lastLabel = unconnectedLabels[-1]
                
                #create an edge from intermediateNode to curNode with unconsumedMatch
                g.add_edge(interNode, curNode, {edgeLabelAttr:(sidx+lcp[idx+1], reverseMatchLen)})
                
                #move to intermediateNode
                curNode = interNode
                
                #change the label of lastNode to lastLabel - unconsumedMatch
                unconnectedLabels[-1] = (lastLabel[0], lastLabel[1]-reverseMatchLen)
        
    return g    
    

g = suffixArrayToSuffixTree(genome1, suffixArray1, lcp1)
#### print edges
##edges = [(e[0], e[1], genome[e[2][edgeLabelAttr][0]:sum(e[2][edgeLabelAttr])]) for e in g.edges_iter(data=True)]
##for e in edges:
##    print(e)
         

## output
with open(outputFile, "w") as f:
    for e in g.edges_iter(data=True):
        #print(genome1[e[2][edgeLabelAttr][0]:sum(e[2][edgeLabelAttr])])
        f.writelines(genome1[e[2][edgeLabelAttr][0]:sum(e[2][edgeLabelAttr])] + '\n')

print('done.')
