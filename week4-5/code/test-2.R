source("code/utility.assembly.R")
source("code/utility.graph.R")

inputFile <- "data/universal_string.txt"
#inputFile <- "C:/Users/ashis/Downloads/dataset_58_10.txt"
edgeFile <- "results/5.3-tmp-edges.txt"
outputFile <- "results/5.3-string-reconstruction-2.txt"

inputs <- readLines(con=inputFile, warn=F)
k <- as.integer(inputs[1])

universalString <- function(k){
  kmers <- sapply(0:(2^k-1), function(n){return(paste(rev(as.integer(intToBits(n))[1:k]), collapse=""))})
  deb <- getDeBrujinGraph(patterns=kmers)
  writeEdgesInFile(from=deb[,"from"], to=deb[,"to"], file=edgeFile, sep=" -> ")
  g <- graphFromFile(inputFile=edgeFile)
  ecycle <- eulerian(g)
  ustr <- substr(genomeFromEulerPath(ecycle), k, 2^k+k-1)
  return(ustr)
}


ustr <- universalString(k)
writeLines(ustr,con=outputFile)