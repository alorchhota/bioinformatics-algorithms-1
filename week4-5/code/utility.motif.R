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

#####################################################
##### convert a string to C-like char vector ########
#####################################################
stringToCharVec <- function(string){
  if(is.null(string) || is.na(string))
    return(c())
  
  charVec <- strsplit(string, split="")[[1]]
  return(charVec)
}

#####################################################
#### calculate hamming distance between 2 strings ###
#####################################################
hammingDistance <- function(string1, string2){
  ch1 <- stringToCharVec(string1)
  ch2 <- stringToCharVec(string2)
  return(sum(ch1!=ch2))
}

#########################################
#### Given a pattern and a dna   ########
#### Find the minimun distance    #######
#### between the pattern and dna    #####
#### when pattern slides over the dna ###
#########################################

distanceBetweenPatternAndDna <- function(pattern, dna){
  k <- str_length(pattern)
  kmers <- getKmers(dna, k=k, unique=T)
  
  distances <- sapply(kmers, hammingDistance, pattern)
  return(min(distances))
}

################################################
#### Given a pattern and a dna collection ######
#### Find the distance between them   ##########
################################################

distanceBetweenPatternAndDnaCollection <- function(pattern, dnaCollection){
  distances <- sapply(dnaCollection, function(dna){
    return(distanceBetweenPatternAndDna(pattern=pattern, dna=dna))
  })
  
  return(sum(distances))
}

#######################################
##### Get all possible DNA ############
#######################################
getAllPossibleDna <- function(len){
  allDnaVector <- permutations(4,len,c("A","T","G","C"), repeats.allowed=T)
  allDna <- apply(allDnaVector, 1, paste0, collapse="")
  return(allDna)
}

###############################################
#### Given a profile matrix           #########
#### Find the probability of a dna segment ####
###############################################
probabilityOfDna <- function(dna, profileMatrix){
  if(is.null(dna) || is.na(dna))
    return(0)
  
  dnaVec <- stringToCharVec(dna)
  prob <- sapply(1:length(dnaVec), function(pos){
    return(profileMatrix[pos, dnaVec[pos]])
  })
  
  probability <- prod(prob)
  return(probability)

}


# Given a profile matrix Profile, we can evaluate the probability of every k-mer in a string Text and find a Profile-most probable k-mer in Text, i.e., a k-mer that was most likely to have been generated by Profile. For example, ACGGGGATTACC is the Profile-most probable 12-mer in GGTACGGGGATTACCT. Indeed, every other 12-mer in this string has probability 0.
# 
# Profile-most Probable k-mer Problem: Find a Profile-most probable k-mer in a string.
# Input: A string Text, an integer k, and a k × 4 matrix Profile.
# Output: A Profile-most probable k-mer in Text.

profileMostProbableKmer <- function(dna, k, profile){
  kmers <- getKmers(text=dna, k=k, unique=T)
  probabilities <- sapply(kmers, probabilityOfDna, profile)
  maxProbableKmerIndex <- which(probabilities==max(probabilities))
  mostProbableKmer <- kmers[maxProbableKmerIndex]
  return(mostProbableKmer)
}

###############################################
#### Given a set of motifs         ############
#### Calculate the profile matrix  ############
###############################################

motifsToProfile_v0.1 <- function(motifs, pseudocount=F){
  motifMatrix <- sapply(motifs, stringToCharVec)
  
  if(!pseudocount){
    profile <- apply(motifMatrix, 1, function(row){
      count <- length(row)
      freqA <- sum(row=="A")/count
      freqC <- sum(row=="C")/count
      freqG <- sum(row=="G")/count
      freqT <- sum(row=="T")/count
      return(c(freqA, freqC, freqG, freqT))
    })
  }
  else{
    profile <- apply(motifMatrix, 1, function(row){
      count <- length(row) + 4
      freqA <- (sum(row=="A")+1)/count
      freqC <- (sum(row=="C")+1)/count
      freqG <- (sum(row=="G")+1)/count
      freqT <- (sum(row=="T")+1)/count
      return(c(freqA, freqC, freqG, freqT))
    })
  }
  
  
  profile <- t(profile)
  colnames(profile) <- c("A","C","G","T")
  return(profile)
}

motifsToProfile <- function(motifs, pseudocount=F){
  noOfMotifs <- length(motifs)
  count <- ifelse(pseudocount, noOfMotifs+4, noOfMotifs)
  motifMatrix <- matrix(unlist(strsplit(x=motifs, split="")), ncol=noOfMotifs)
  
  if(!pseudocount){
    profile <- apply(motifMatrix, 1, function(row){
      freqA <- sum(row=="A")/count
      freqC <- sum(row=="C")/count
      freqG <- sum(row=="G")/count
      freqT <- sum(row=="T")/count
      return(c(freqA, freqC, freqG, freqT))
    })
  }
  else{
    profile <- apply(motifMatrix, 1, function(row){
      freqA <- (sum(row=="A")+1)/count
      freqC <- (sum(row=="C")+1)/count
      freqG <- (sum(row=="G")+1)/count
      freqT <- (sum(row=="T")+1)/count
      return(c(freqA, freqC, freqG, freqT))
    })
  }
  
  
  profile <- t(profile)
  colnames(profile) <- c("A","C","G","T")
  return(profile)
}

###############################################
#### Given a set of motifs         ############
#### Calculate the count matrix    ########
###############################################

motifsToCountMatrix <- function(motifs){
  noOfMotifs <- length(motifs)
  motifMatrix <- matrix(unlist(strsplit(x=motifs, split="")), ncol=noOfMotifs)
  
  countMat <- apply(motifMatrix, 1, function(row){
    freqA <- sum(row=="A")
    freqC <- sum(row=="C")
    freqG <- sum(row=="G")
    freqT <- sum(row=="T")
    return(c(freqA, freqC, freqG, freqT))
  })
  
  countMat <- t(countMat)
  
  colnames(countMat) <- c("A","C","G","T")
  return(countMat)
}

countMatrixToProfile <- function(countMatrix, pseudocount=F){
  if(pseudocount){
    countMatrix <- countMatrix + 1
  }
  
  noOfMotifs <- sum(countMatrix[1,])
  profile <- countMatrix/noOfMotifs
  return(profile)
}

updateCountMatrixForMotifDeletion <- function(countMatrix, motif){
  motifCount <- motifsToCountMatrix(motif)
  updatedCountMatrix <- countMatrix - motifCount
  return(updatedCountMatrix)
}

updateCountMatrixForMotifAddition <- function(countMatrix, motif){
  motifCount <- motifsToCountMatrix(motif)
  updatedCountMatrix <- countMatrix + motifCount
  return(updatedCountMatrix)
}

###############################################
#### Calculate score of motifs  ###############
###############################################
scoreOfMotifs <- function(motifs){
  profile <- motifsToProfile(motifs)
  noOfMotifs <- length(motifs)
  score <- scoreOfMotifsFromProfile(profile, noOfMotifs)
  return(score)
}

scoreOfMotifsFromProfile <- function(profile, noOfMotifs){
  posScores <- apply(profile, 1, function(row){
    frequentColIndex <- which(row==max(row))[1]
    onePosScore <- sum(row[-frequentColIndex])
    return(onePosScore)
  })
  
  score <- sum(posScores) * noOfMotifs
  return(score)
}

scoreOfMotifsFromCountMatrix <- function(countMatrix){
  posScores <- apply(countMatrix, 1, function(row){
    frequentColIndex <- which(row==max(row))[1]
    onePosScore <- sum(row[-frequentColIndex])
    return(onePosScore)
  })
  
  score <- sum(posScores)
  return(score)
}


###############################################
#### Find Consensus  motif  ###################
###############################################
consensusMotif <- function(motifs){
  profile <- motifsToProfile(motifs)
  consensusMotif <- consensusMotifFromProfile(profile)
  return(consensusMotif)
}

consensusMotifFromProfile <- function(profile){
  cols <- colnames(profile)
  colIndexes <- apply(profile, 1, function(row){
    return(which(row==max(row))[1])
  })
  consensusVec <- cols[colIndexes]
  consensus <- paste0(consensusVec, collapse="")
  return(consensus)
}

###############################################
#### Calculate motifs from profile   ##########
#### motif = Profile-most probable kmer #######
###############################################
profileToMotifs <- function(profile, dnaStrings){
  k <- nrow(profile)
  motifs <- sapply(dnaStrings, function(dna){
    motif <- profileMostProbableKmer(dna=dna, k=k, profile=profile)[1]
    return(motif)
  })
  
  return(motifs)
}


profileRandomKmer <- function(profile, kmers, dna=NA, k=0){
  ## if kmers are not gives, then kmers are generated from dna and k
  if(is.null(kmers) || is.na(kmers))
    kmers <- getKmers(text=dna, k=k, unique=T)
  
  probabilities <- sapply(kmers, probabilityOfDna, profile)
  randomKmer <- sample(kmers, size=1, prob=probabilities)
  return(randomKmer)
}

#debug(probabilityOfDna)


composition <- function(text, k){
  kmers <- getKmers(text=text, k=k, unique=F)
  sortedKmers <- sort(kmers)
  return(sortedKmers)
}


getOverlapGraph <- function(patterns){
  
  patterns <- sort(patterns)
  names(patterns) <- patterns
  
  edges <- c()
  addAdjancency <- function(pattern){
    lastPart <- substr(pattern,start=2,stop=str_length(pattern))
    candidatePatterns <- paste0(lastPart, c("A","C","G","T"))
    candidateValues <- patterns[candidatePatterns]
    neighbors <- candidateValues[!is.na(candidateValues)]
    
    lapply(neighbors, function(neighbor){
      edges <<- append(edges, paste0(pattern, " -> ", neighbor))
      return(NA)
    })
    
    return(NA)
  }
  
  sapply(patterns, addAdjancency)
  
  return(edges)
}


#debug(getOverlapGraph)


