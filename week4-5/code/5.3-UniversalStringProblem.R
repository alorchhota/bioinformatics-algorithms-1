source("code/utility.assembly.R")
source("code/utility.graph.R")

inputFile <- "data/universal_string.txt"
#inputFile <- "C:/Users/ashis/Downloads/dataset_58_10.txt"
edgeFile <- "results/5.3-tmp-edges.txt"
outputFile <- "results/5.3-string-reconstruction.txt"

inputs <- readLines(con=inputFile, warn=F)
k <- as.integer(inputs[1])

universalString <- function(k){
  #t1 <- Sys.time()
  #print("creating kmers")
  kmers <- sapply(0:(2^k-1), function(n){return(paste(rev(as.integer(intToBits(n))[1:k]), collapse=""))})
  #print("creating debrujin graph")
  deb <- getDeBrujinGraph(patterns=kmers)
  #print("making graph")
  writeEdgesInFile(from=deb[,from], to=deb[,to], file=edgeFile, sep=" -> ")
  g <- graphFromFile(inputFile=edgeFile)
  #t2 <- Sys.time()
  #print(t2-t1)
  #print("euler cycle")
  start <- paste0(rep("0", k-1), collapse="")
  ecycle <- eulerian(g, start)
  #print("creating universal string")
  ustr <- substr(genomeFromEulerPath(ecycle), k, 2^k+k-1)
  return(ustr)
}

t1 = Sys.time()
ustr <- universalString(k)
print("write string in file")
writeLines(ustr,con=outputFile)
t2 = Sys.time()
print(t2-t1)


# k <- 19
# t1 <- Sys.time()
# kmers <- sapply(0:(2^k-1), function(n){return(paste(rev(as.integer(intToBits(n))[1:k]), collapse=""))})
# deb <- getDeBrujinGraph(patterns=kmers)
# t2 <- Sys.time()
# print(t2-t1)