# Frequent Words with Mismatches and Reverse Complements: Find the most frequent k-mers in a string.
# Input: A string Text and an integer k.
# Output: All most frequent k-mers in Text.

library("stringr")
library("gtools")

strReverse <- function(x){sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")}

complement <- function(str){
  strVec <- strsplit(str, NULL)[[1]]
  complementVec <- rep('',length(strVec))
  complementVec[which(strVec=="A")] <- "T"
  complementVec[which(strVec=="T")] <- "A"
  complementVec[which(strVec=="G")] <- "C"
  complementVec[which(strVec=="C")] <- "G"
  return(paste0(complementVec, collapse=""))
}

reverseComplement <- function (pattern){
  revPattern <- strReverse(pattern)  
  revCom <- complement(revPattern)
  return(revCom)
}


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
  
  sapply(d:d, generateApproxKmersWithNChanges)
  return(unique(approxKmers))
}


freuent.k.mers.with.mismatch <- function(text, k, d){
  
  kmers <- sapply(1:(str_length(text)-k+1), function(start){
    end <- start + k -1
    si <- substr(x=text, start=start, stop=end)
    return(si)
  })
  
  freq <- c()
  allApproxKmers <- c()
  
  # foreach kmer, save all possible approximate kmers 
  findAndSaveApproximateKmers <- function(kmer){
	approxKmers <- getApproximateKmers(kmer, d)
    allApproxKmers <<- c(allApproxKmers, approxKmers)
    return(NA)
  }
  sapply(kmers, findAndSaveApproximateKmers)
  
  # foreach reverse complement, save all possible approximate kmers
  rckmers <- sapply(kmers, reverseComplement)
  sapply(rckmers, findAndSaveApproximateKmers)
  
  # calculate frequency and find max
  freq <- table(allApproxKmers)
  maxFreq <- max(freq)
  maxFreqItems <- names(which(freq==maxFreq))
  return(list(kmers=maxFreqItems, max.frequency=maxFreq, allFreq=freq))
  
}

text <- "ACGTTGCATGTCGCATGATGCATGAGAGCT"
k <- 4
d <- 1
frequentKmers <- freuent.k.mers.with.mismatch(text=text, k=k, d=d)
print(frequentKmers$kmers)

