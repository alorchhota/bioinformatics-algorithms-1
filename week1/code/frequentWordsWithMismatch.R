# Frequent Words Problem: Find the most frequent k-mers in a string.
# Input: A string Text and an integer k.
# Output: All most frequent k-mers in Text.

library("stringr")
library("gtools")

getApproximateKmers <- function(kmer, d){
  kmerVec <- strsplit(kmer, "")[[1]]
  approxKmers <- c()
  
  generateApproxKmersWithNChanges <- function(nChangedPos){
    changedPos <- combn(1:length(kmerVec), nChangedPos)
    changedLetters <- permutations(4,nChangedPos,c("A","T","G","C"), repeats.allowed=T)
    
    sapply(1:ncol(changedPos), function(changePosNo){
      sapply(1:nrow(changedLetters), function(changeLettersNo){
        approxKmerVec <- kmerVec
        approxKmerVec[changedPos[,changePosNo]] <- changedLetters[changeLettersNo,]
        approxKmer <- paste0(approxKmerVec, collapse="")
        if(is.null(approxKmers) ||  is.na(approxKmers[approxKmer]))
          approxKmers[approxKmer] <<- approxKmer
        return(NA)
      })
      
      return(NA)
    })
    
    
    return(NA)
  }
  
  #debug(generateApproxKmersWithNChanges)
  
  sapply(d:d, generateApproxKmersWithNChanges)
  return(unique(approxKmers))
  
}

#debug(getApproximateKmers)

freuent.k.mers.with.mismatch <- function(text, k, d){
  
  kmers <- sapply(1:(str_length(text)-k+1), function(start){
    end <- start + k -1
    si <- substr(x=text, start=start, stop=end)
    return(si)
  })
  
  freq <- c()
  allApproxKmers <- c()
  
  # foreach kmer update frequency of possible approximate kmers 
  
  countKmer <- 0
  
  updateFrequencyOfApproximateKmers <- function(kmer){
    countKmer <<- countKmer + 1
    print(countKmer)
    
    approxKmers <- getApproximateKmers(kmer, d)
    #approxKmers <- unique(approxKmers)
    allApproxKmers <<- c(allApproxKmers, approxKmers)
    
#     updateFrequencyOfApproximateKmer <- function(approxKmer){
#       if(is.null(freq) || is.na(freq[approxKmer])){
#         freq[approxKmer] <<- 1
#       }
#       else{
#         freq[approxKmer] <<- freq[approxKmer] + 1
#       }
#       return(NA)
#     }
    
    
    
#    sapply(approxKmers, updateFrequencyOfApproximateKmer)
    return(NA)
  }
  
  
  print(paste("No Of Kmers:", length(kmers)))
  sapply(kmers, updateFrequencyOfApproximateKmers)
  
    
  # find the frequent items
#   maxFreq <- max(freq)
#   maxFreqItems <- names(freq[freq==maxFreq])
#   return(list(kmers=maxFreqItems, max.frequency=maxFreq))
  
  freq <- table(allApproxKmers)
  maxFreq <- max(freq)
  maxFreqItems <- names(which(freq==maxFreq))
  return(list(kmers=maxFreqItems, max.frequency=maxFreq, allFreq=freq))
  
}
#debug(freuent.k.mers.with.mismatch)


#debug(freuent.k.mers.with.mismatch)

# text <- "ACGTTGCATGTCGCATGATGCATGAGAGCT"
# k <- 4
# d <- 1
# frequentKmers <- freuent.k.mers.with.mismatch(text=text, k=k, d=d)
# print(frequentKmers$kmers)

# text <- "CACAGTAGGCGCCGGCACACACAGCCCCGGGCCCCGGGCCGCCCCGGGCCGGCGGCCGCCGGCGCCGGCACACCGGCACAGCCGTACCGGCACAGTAGTACCGGCCGGCCGGCACACCGGCACACCGGGTACACACCGGGGCGCACACACAGGCGGGCGCCGGGCCCCGGGCCGTACCGGGCCGCCGGCGGCCCACAGGCGCCGGCACAGTACCGGCACACACAGTAGCCCACACACAGGCGGGCGGTAGCCGGCGCACACACACACAGTAGGCGCACAGCCGCCCACACACACCGGCCGGCCGGCACAGGCGGGCGGGCGCACACACACCGGCACAGTAGTAGGCGGCCGGCGCACAGCC"
# k <- 10
# d <- 2
# frequentKmers <- freuent.k.mers.with.mismatch(text=text, k=k, d=d)
# print(frequentKmers$kmers)


inputFile <- "data/frequent_words_mismatch_data4.txt"
inputs <- readLines(con=inputFile, n=1)
inputs <- strsplit(inputs, " ")[[1]]
text <- inputs[1]
k <- as.integer(inputs[2])
d <- as.integer(inputs[3])
frequentKmers <- freuent.k.mers.with.mismatch(text=text, k=k, d=d)
cat(frequentKmers$kmers)


# app <- getApproximateKmers("AT",2)
# print(length(app))
# print(app)