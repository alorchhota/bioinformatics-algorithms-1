source("code/utility.motif.R")

#inputFile <- "data/overlap_graph_1-1.txt"
inputFile <- "C:\\Users\\Ashis\\Downloads\\dataset_53_7.txt"
outputFile <- "results/3.2-overlap-graph.txt"


inputs <- readLines(con=inputFile, warn=F)
patterns <- inputs


print(Sys.time())
overlapGraph <- getOverlapGraph(patterns=patterns)
write.table(x=overlapGraph, file=outputFile, quote=F, sep="\n", row.names=F, col.names=F)