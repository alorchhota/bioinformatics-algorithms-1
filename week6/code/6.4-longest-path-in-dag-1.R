source("code/utility.graph.R")
inputFile <- "data/6.4-longest-path-in-dag-1.txt"
#inputFile <- "C:/Users/Ashis/Downloads/dataset_75_5 (4).txt"
outputFile <- "results/6.4-longest-path-in-dag-1.txt"

inputs <- readLines(con=inputFile, warn=F)
start <- inputs[1]
end <- inputs[2]
edges <- inputs[3:length(inputs)]

g <- weightedGraphFromEdges(edges)

getLongestPathInDAG(g, "0","4")
