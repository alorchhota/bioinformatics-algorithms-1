source("code/utility.assembly.R")
source("code/utility.graph.R")

inputFile <- "data/str_recon-2.txt"
#inputFile <- "C:/Users/ashis/Downloads/dataset_58_6.txt"
outputFile <- "results/5.2-string-reconstruction.txt"

gR <- graphFromFile(inputFile)
print(numEdges(gR))
epath <- eulerian(gR)
genome <- genomeFromEulerPath(epath)
writeLines(genome,con=outputFile)