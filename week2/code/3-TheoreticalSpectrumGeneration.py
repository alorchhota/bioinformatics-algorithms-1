import os
import csv
import sys
import re
import importlib

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
inputFile = './data/3.theoretical-spectrum-data-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_98_3.txt'
outputFile = './results/3.theoretical-spectrum.txt'
aminoAcidMassMapFile = './data/integer_mass_table.txt'


# set current directory
os.chdir(curDir)

# read input
with open(inputFile) as f:
    peptide = f.readline().strip()

# load SeqUtil methods
sys.path.append('code')
import SeqUtil
importlib.reload(SeqUtil)

### original work ####
# load amion acid mass
AAMassMap = SeqUtil.loadAminoAcidMass(aminoAcidMassMapFile)

# generate spectrum
spectrum = SeqUtil.cyclospectrum(peptide, AAMassMap)

# output
with open(outputFile, "w") as f:
    f.writelines(" ".join([str(s) for s in spectrum]))


