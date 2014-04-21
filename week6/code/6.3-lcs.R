source("code/utility.string.R")
inputFile <- "data/6.3-lcs-1.txt"
#inputFile <- "C:/Users/Ashis/Downloads/dataset_75_5 (4).txt"
outputFile <- "results/6.3-lcs.txt"

inputs <- readLines(con=inputFile, warn=F)
s <- inputs[1]
t <- inputs[2]

seq <- lcs(s,t)
cat(seq,file=outputFile)
