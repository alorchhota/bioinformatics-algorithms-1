# Clump Finding Problem: Find patterns forming clumps in a string.
# Input: A string Genome, and integers k, L, and t.
# Output: All distinct k-mers forming (L, t)-clumps in Genome.

library("stringr")

getClumps <- function(genome, k, L, t){
  lenGenome <- str_length(genome)
  
  findAllKmers <- function(string, k, unique=F){
    len <- str_length(string)
    
    kmers <- sapply(1:(len-k+1), function(i){
      return(substr(string, i, i+k-1))
    })
    
    if(unique)
      return(unique(kmers))
    else
      return(kmers)
    
  }
  
  findFrequentItems <- function(items, minFreq=1){
    freq <- table(items)
    freqItems <- names(freq[freq>=minFreq])
    return(freqItems)
  }
  
  
  clumps <- c()
  saveClump <- function(item){
    if(!any(clumps==item))
      clumps <<- c(clumps, item)
  }
  
  tmp <- sapply(1:(lenGenome-L+1), function(start){
    end <- start + L - 1
    windowGenome <- substr(genome, start, end)
    
    # find all k-mers in a window
    kmersInWindow <- findAllKmers(windowGenome, k=k, unique=F)
    freqItems <- findFrequentItems(kmersInWindow, minFreq=t)
    if(any(freqItems=="GAAGA")) 
      debug(updateClumpFreq)
    sapply(freqItems, saveClump)
    return(NA)
  })
  
  return(clumps)
  
}

inputFile <- "data/ClumpFinding3.txt"
inputs <- readLines(con=inputFile, n=2)
genome <- inputs[1]
numbers <- as.integer(strsplit(inputs[2], " ")[[1]])
k <- numbers[1]
L <- numbers[2]
t <- numbers[3]

# genome <- "CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA"
# k <- 5
# L <- 50
# t <- 4

clumps <- getClumps(genome, k, L, t)
cat(clumps)
