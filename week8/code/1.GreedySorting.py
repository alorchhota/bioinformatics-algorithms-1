import os
import csv
import sys
import re
import importlib

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
inputFile = './data/1.greedy_sorting-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_286_3.txt'
outputFile = './results/1.greedy_sorting-1.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

perm = inputs[0].strip()[1:-1]
perm = [int(x) for x in re.split('\\s+', perm)]
print(perm)

n = len(perm)
with open(outputFile, "w") as f:
    for i in range(n):
        if perm[i] != i+1 or perm[i] == - (i+1):
            if i+1 in perm:
                j = perm.index(i+1)
            else:
                j = perm.index(-i-1)

            # reverse with sign inversing
            perm[i:j+1] = [-x for x in perm[j:i:-1]+[perm[i]]]
            f.writelines('(' + ' '.join([ ['-','+'][int(x>0)] + str(abs(x)) for x in perm]) + ')\n')
            
        if perm[i] != i+1:
            # reverse only the first one
            perm[i] = - perm[i]
            f.writelines('(' + ' '.join([ ['-','+'][int(x>0)] + str(abs(x)) for x in perm]) + ')\n')
        

print(perm)


# output
#
#    
