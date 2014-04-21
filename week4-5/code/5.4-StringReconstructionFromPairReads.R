source("code/utility.assembly.R")
source("code/utility.graph.R")

#inputFile <- "data/pair_end-2.txt"
inputFile <- "C:/Users/ashis/Downloads/dataset_59_14.txt"
outputFile <- "results/5.4-string-reconstruction.txt"
edgeFile <- "results/5.4-tmp-edges.txt"

inputs <- readLines(con=inputFile,warn=F)
d <- as.integer(inputs[1])
pairReads <- inputs[-1]


deb <- getDeBrujinGraphFromPairedReads(pairReads, sep="-")
writeEdgesInFile(from=deb[,from], to=deb[,to], file=edgeFile, sep=" -> ")
gR <- graphFromFile(edgeFile)
#print(numEdges(gR))
epath <- eulerian(gR)
genome <- genomeFromEulerPathOfPairReads(epath, d)
writeLines(genome,con=outputFile)