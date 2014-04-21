import os
import csv
import sys
import re
import importlib

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
#curDir = 'D:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
inputFile = './data/2.peptide_encoding_data-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_96_6.txt'
outputFile = './results/2.peptide_encoding.txt'
codonMapFile = './data/RNA_codon_table_1.txt'

# set current directory
os.chdir(curDir)

# read input
with open(inputFile) as f:
    dna = f.readline().strip()
    target = f.readline().strip()

# load SeqUtil methods
sys.path.append('code')
import SeqUtil
importlib.reload(SeqUtil)


lenTargetPeptide = len(target)
lenTargetDna = lenTargetPeptide*3

foundDnas = []  # storage of target dan sequences

## search amino acids in forward direction
for start in range(0,3):                # loop for reading frame start
    curDna = dna[start:]                
    curRna = SeqUtil.dna2rna(curDna)
    peptides = SeqUtil.rna2peptide(curRna, mapFile='./data/RNA_codon_table_1.txt')

    pepStartIndex = 0
    for pep in peptides:
        # find positions of target in peptide
        targetPositions = [pos for pos in range(0, len(pep)) if pep[pos:(pos+lenTargetPeptide)]==target]

        # get corresponding dna sequence and save
        targetDnas = [curDna[(pepStartIndex+pos*3):(pepStartIndex+pos*3+lenTargetDna)] for pos in targetPositions]
        foundDnas.extend(targetDnas)

        # update pepStartIndex for the next peptide
        pepStartIndex = pepStartIndex + (len(pep)+1)*3


## search amino acids in reverse direction
complDna = SeqUtil.reverseComplement(dna)
for start in range(0,3):                # loop for reading frame start
    curDna = complDna[start:]
    curRna = SeqUtil.dna2rna(curDna)
    peptides = SeqUtil.rna2peptide(curRna, mapFile='./data/RNA_codon_table_1.txt')

    pepStartIndex = 0
    for pep in peptides:
        # find positions of target in peptide
        targetPositions = [pos for pos in range(0, len(pep)) if pep[pos:(pos+lenTargetPeptide)]==target]

        # get corresponding dna sequence and save
        targetDnas = [curDna[(pepStartIndex+pos*3):(pepStartIndex+pos*3+lenTargetDna)] for pos in targetPositions]
        targetDnas = [SeqUtil.reverseComplement(dna) for dna in targetDnas]    # reverse complement
        foundDnas.extend(targetDnas)

        # update pepStartIndex for the next peptide
        pepStartIndex = pepStartIndex + (len(pep)+1)*3


# output
with open(outputFile, "w") as f:
    f.writelines("\n".join(foundDnas))
