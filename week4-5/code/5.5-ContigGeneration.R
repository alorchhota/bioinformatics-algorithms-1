source("code/utility.assembly.R")

inputFile <- "data/contig_generation-2.txt"
inputFile <- "C:/Users/ashis/Downloads/dataset_60_5.txt"
outputFile <- "results/5.5-contig-generation.txt"
edgeFile <- "results/5.5-tmp-edges.txt"

inputs <- readLines(con=inputFile,warn=F)
patterns <- inputs

deb <- getDeBrujinGraph(patterns=patterns)
writeEdgesInFile(from=deb[,from], to=deb[,to], file=edgeFile, sep=" -> ")
g <- graphFromFile(inputFile=edgeFile)
cons <- contigs(g)
writeLines(cons,con=outputFile,sep="\n")

