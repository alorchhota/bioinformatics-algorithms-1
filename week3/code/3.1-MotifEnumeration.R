inputFile <- "data/motif_enumeration_data-4.txt"
inputs <- readLines(con=inputFile, warn=F)
numbers <- strsplit(inputs[1], " ")[[1]]
k <- as.integer(numbers[1])
d <- as.integer(numbers[2])
dnaCollection <- inputs[2:length(inputs)]

library("stringr")
library("gtools")

###########################################################################
#### Given a kmer, Find all kmers differing by at most 'd' mutations ######
###########################################################################
getApproximateKmers <- function(kmer, d){
  kmerVec <- strsplit(kmer, "")[[1]]
  
  changedPos <- combn(1:length(kmerVec), d)
  changedLetters <- permutations(4,d,c("A","T","G","C"), repeats.allowed=T)
  
  allPositionsApproxKmers <- sapply(1:ncol(changedPos), function(changePosNo){
    onePositionApproxKmers <- sapply(1:nrow(changedLetters), function(changeLettersNo){
      approxKmerVec <- kmerVec
      approxKmerVec[changedPos[,changePosNo]] <- changedLetters[changeLettersNo,]
      approxKmer <- paste0(approxKmerVec, collapse="")
      return(approxKmer)
    })
    
    return(onePositionApproxKmers)
  })
  
  return(unique(as.vector(allPositionsApproxKmers)))
}

###########################################################################
###### given a dna (text), Find all possible kmers ########################
###########################################################################
getKmers <- function(text, k, unique=F){
  kmers <- sapply(1:(str_length(text)-k+1), function(start){
    end <- start + k -1
    si <- substr(x=text, start=start, stop=end)
    return(si)
  })
  
  if(unique)
    return(unique(kmers))
  else
    return(kmers)
}

##################################################################################
###### does two string approximately match (differ by at most d positions) #######
##################################################################################
match <- function(string1, string2, d=0){
  charStr1 <- strsplit(x=string1, split="")[[1]]
  charStr2 <- strsplit(x=string2, split="")[[1]]
  diff <- sum( charStr1 != charStr2)
  return(diff <= d)
}


#####################################################
###### Given a dna (string), a segment (pattern) ####
###### Find all the positions in the dna   ##########
###### where the segment matches the dna   ##########
###### (at most d mutations is acceptable) ##########
#####################################################
findMatches <- function(string, pattern, d=0, indexBase=1){
  
  lenPattern <- str_length(pattern)
  lenStr <- str_length(string)
  
  if(lenPattern > lenStr)
    return(c())
  
  charPattern <- strsplit(x=pattern, split="")[[1]]
  
  matches <- sapply(1:(lenStr-lenPattern+1), function(start){
    end <- start + lenPattern - 1
    curStr <- substr(string, start, end)
    charCurStr <- strsplit(x=curStr, split="")[[1]]
    diff <- sum( charPattern != charCurStr)
    return(diff <= d)
  })
  
  matchIndexes <- which(matches==T)
  if(indexBase != 1)
    matchIndexes <- matchIndexes - 1 + indexBase
  
  return(matchIndexes)
}




# MOTIFENUMERATION(Dna, k, d)
# for each k-mer a in Dna
# for each k-mer a' differing from a by at most d mutations
# if a' appears in each string from Dna with at most d mutations
# output a'

motifEnumerate <- function(dnaCollection, k, d){
  
  allKmers <- lapply(dnaCollection, getKmers, k)
  allKmers <- unique(unlist(allKmers))
  
  allApproximateKmers <- lapply(allKmers, getApproximateKmers, d)
  allApproximateKmers <- unique(unlist(allApproximateKmers))
  
  hasSegmentInDna <- function(dna, segment, d=0){
    matchIndexes <- findMatches(string=dna, pattern=segment, d=d)
    return(length(matchIndexes) > 0)
  }
  
  hasSegmentInAllDna <- function(dnaCollection, segment, d=0){
    hasInDnaCollection <- sapply(dnaCollection, hasSegmentInDna, segment, d)
    return(all(hasInDnaCollection))
  }
  
  isMotif <- sapply(allApproximateKmers, function(kmer){
    return(hasSegmentInAllDna(dnaCollection, kmer, d))
  })
  
  motifs <- allApproximateKmers[isMotif]
  
  return(motifs)
  
}

#debug(motifEnumerate)


motifs <- motifEnumerate(dnaCollection, k, d)
cat(motifs)