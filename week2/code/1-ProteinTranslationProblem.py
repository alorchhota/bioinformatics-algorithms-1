import os
import csv

# settings
curDir = 'E:/Copy/Coursera/Bioinformatics Algorithms (part-I)/MyPrograms/week2'
inputFile = './data/1.protein_translate_data-1.txt'
#inputFile = 'C:/Users/Ashis/Downloads/dataset_96_3.txt'
outputFile = './results/1.protein_translate.txt'
codonMapFile = './data/RNA_codon_table_1.txt'

# set current directory
os.chdir(curDir)

# read codon map file and create a codon map
with open(codonMapFile) as f:
    freader = csv.reader(f, delimiter=' ')
    codonMap = {row[0]:row[1] for row in freader}

# read input and translate to peptide
with open(inputFile) as f:
	pattern = f.readline()
codons = [pattern[start:(start+3)] for start in range(0, len(pattern)-3, 3)]
peptide = "".join([codonMap[c] for c in codons])

# output
with open(outputFile, "w") as f:
    f.write(peptide)
