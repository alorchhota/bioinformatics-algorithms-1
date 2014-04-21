import os
import csv
import sys
import re
import importlib

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
inputFile = './data/7.convolutional_seq_data-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_104_7.txt'
outputFile = './results/7.convolutional_seq.txt'


# set current directory
os.chdir(curDir)

# read input
with open(inputFile) as f:
    M = int(f.readline().strip())
    N = int(f.readline().strip())
    spectrumLine = f.readline().strip()
    spectrum = [int(n) for n in spectrumLine.split(" ")]


# load SeqUtil methods
sys.path.append('code')
import SeqUtil
importlib.reload(SeqUtil)

### original work ####
# generate sequence from spectrum
seq = SeqUtil.convolutionCyclopeptideSequencing(spectrum, M, N)

# output
with open(outputFile, "w") as f:
    f.writelines("-".join([str(s) for s in seq]))


