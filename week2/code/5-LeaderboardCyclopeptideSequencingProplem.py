import os
import csv
import sys
import re
import importlib

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
inputFile = './data/5.leaderboard_data-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_102_4 (1).txt'
outputFile = './results/5.leaderboard.txt'
aminoAcidMassMapFile = './data/integer_mass_table.txt'


# set current directory
os.chdir(curDir)

# read input
with open(inputFile) as f:
    N = int(f.readline().strip())
    spectrumLine = f.readline().strip()
    spectrum = [int(n) for n in spectrumLine.split(" ")]


# load SeqUtil methods
sys.path.append('code')
import SeqUtil
importlib.reload(SeqUtil)

### original work ####
# load amion acid mass
AAMassMap = SeqUtil.loadAminoAcidMass(aminoAcidMassMapFile)

# generate sequence from spectrum
seq = SeqUtil.leaderboardCyclopeptideSequencing(spectrum, N, AAMassMap.values())

# output
with open(outputFile, "w") as f:
    f.writelines("-".join([str(s) for s in seq]))


