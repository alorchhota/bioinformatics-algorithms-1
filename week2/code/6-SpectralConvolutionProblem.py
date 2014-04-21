import os
import csv
import sys
import re
import importlib

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
inputFile = './data/6.spectral_convolution_data-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_104_4.txt'
outputFile = './results/6.spectral_convolution.txt'
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
conv = SeqUtil.spectralConvolution(spectrum)
conv = sorted(conv)

# output
with open(outputFile, "w") as f:
    f.writelines(" ".join([str(m) for m in conv]))


