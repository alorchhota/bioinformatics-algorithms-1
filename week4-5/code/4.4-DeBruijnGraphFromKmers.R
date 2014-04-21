source("code/utility.assembly.R")

inputFile <- "data/debruijn_graph_kmers-1.txt"
#inputFile <- "C:/Users/ashis/Downloads/dataset_55_7 (3).txt"
outputFile <- "results/4.4-debrujin-graph.txt"


inputs <- readLines(con=inputFile, warn=F)
kmers <- inputs
debrujinEdges <- getDeBrujinGraph(kmers)
writeEdgesInFile(from=debrujinEdges[,"from",drop=F], to=debrujinEdges[,"to",drop=F], sep=" -> ", file=outputFile)

