source("code/utility.motif.R")

#inputFile <- "data/string_com-2.txt"
inputFile <- "C:/Users/Ashis/Downloads/dataset_51_3 (1).txt"
outputFile <- "results/3.1-composition-output.txt"


inputs <- readLines(con=inputFile, warn=F)
numbers <- strsplit(inputs[1], " ")[[1]]
k <- as.integer(numbers[1])
text <- inputs[2]

compos <- composition(text=text, k=k)
cat(compos, sep="\n", file=outputFile)