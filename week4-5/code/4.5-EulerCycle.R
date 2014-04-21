source("code/utility.assembly.R")
source("code/utility.graph.R")

library(PairViz)

#inputFile <- "data/eulerian_cycle-2.txt"
inputFile <- "C:/Users/ashis/Downloads/dataset_58_2.txt"
outputFile <- "results/4.5-euler-cycle.txt"

gR <- graphFromFile(inputFile)
print(numEdges(gR))
print(Sys.time())
eulerCycle <- eulerCycle(gR)
print(Sys.time())
writeLines(paste(eulerCycle, collapse="->"),con=outputFile)