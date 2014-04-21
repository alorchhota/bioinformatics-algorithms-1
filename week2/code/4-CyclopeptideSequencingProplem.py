import os
import csv
import sys
import re
import importlib

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
inputFile = './data/4.cycloseq_data-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_100_4.txt'
outputFile = './results/4.cycloseq.txt'
aminoAcidMassMapFile = './data/integer_mass_table.txt'


# set current directory
os.chdir(curDir)

# read input
with open(inputFile) as f:
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
sequences = SeqUtil.cyclopeptideSequencing(spectrum, AAMassMap.values())

# output
with open(outputFile, "w") as f:
    seqStrings = []
    for seq in sequences:
        seqStrings.append("-".join(str(mass) for mass in seq))
    f.writelines(" ".join([s for s in seqStrings]))


