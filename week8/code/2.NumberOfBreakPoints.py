import os
import csv
import sys
import re
import importlib

# settings
#curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week8'
inputFile = './data/2.NumberOfBreakPoints-2.txt'
inputFile = 'C:/Users/Ashis/Downloads/dataset_287_4.txt'
outputFile = './results/2.NumberOfBreakPoints.txt'

# set current directory
os.chdir(curDir)

## read input
with open(inputFile) as f:
    inputs = f.readlines()

perm = inputs[0].strip()[1:-1]
perm = [int(x) for x in re.split('\\s+', perm)]

## add leading 0 and trailing length
perm.insert(0,0)
perm.append(len(perm))
#print(perm)

breaks = [i for i in range(1,len(perm)) if perm[i-1]!=perm[i]-1]
nbreaks = len(breaks)
print(nbreaks)


# output
with open(outputFile, 'w') as f:
    f.writelines(str(nbreaks)+"\n")
