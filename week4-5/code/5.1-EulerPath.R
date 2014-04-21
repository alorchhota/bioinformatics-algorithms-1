source("code/utility.assembly.R")
source("code/utility.graph.R")

#inputFile <- "data/eulerian_path-1.txt"
inputFile <- "C:/Users/ashis/Downloads/dataset_58_5.txt"
outputFile <- "results/5.1-euler-path.txt"

gR <- graphFromFile(inputFile)
print(numEdges(gR))
print(Sys.time())
eulerCycle <- eulerian(gR)
print(Sys.time())
writeLines(paste(eulerCycle, collapse="->"),con=outputFile)