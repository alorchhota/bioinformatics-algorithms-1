source("code/utility.assembly.R")

#inputFile <- "data/debruijn_graph_string-2.txt"
inputFile <- "C:/Users/Ashis/Downloads/dataset_54_6.txt"
outputFile <- "results/4.3-debrujin-graph.txt"


inputs <- readLines(con=inputFile, warn=F)
k <- as.integer(inputs[1])
text <- inputs[2]

comp <- composition(text, k)
debrujinEdges <- getDeBrujinGraph(comp)

writeEdgesInFile(from=debrujinEdges[,"from",drop=F], to=debrujinEdges[,"to",drop=F], sep=" -> ", file=outputFile)


