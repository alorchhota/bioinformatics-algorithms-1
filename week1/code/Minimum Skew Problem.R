# Minimum Skew Problem: Find a position in a genome minimizing the skew.
# Input: A DNA string Genome.
# Output: All integer(s) i minimizing Skew(Prefixi (Text)) among all values of i (from 0 to |Genome|).

library("stringr")

minSkew <- function(genome){
  curSkew = 0
  
  lenGenome <- str_length(genome)
  isDecreasing <- F
  genome <- strsplit(genome, NULL)[[1]]
  
  skews <- sapply(genome, function(letter){
    if(letter=="G"){
      curSkew <<- curSkew + 1
    }
    else if(letter=="C"){
      curSkew <<- curSkew - 1
    }
    return(curSkew)
  })
  
  minSkew = min(skews)
  minSkewIndexes <- which(skews==minSkew, arr.ind=T)
  return(list(minSkew=minSkew, indexes=minSkewIndexes))
}


inputFile <- "data/minimum_skew3.txt"
inputs <- readLines(con=inputFile, n=1)
genome <- inputs[1]

#genome <- "TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT"
minSkew <- minSkew(genome);
cat(minSkew$indexes)