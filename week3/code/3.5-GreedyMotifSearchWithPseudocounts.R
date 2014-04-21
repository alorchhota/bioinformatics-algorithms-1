source("code/utility.motif.R")

inputFile <- "data/greedy_pseudo-6.txt"
inputs <- readLines(con=inputFile, warn=F)
numbers <- strsplit(inputs[1], " ")[[1]]
k <- as.integer(numbers[1])
t <- as.integer(numbers[2])
dnaStrings <- inputs[2:length(inputs)]


# profile <- motifsToProfile(dnaStrings)
# score <- scoreOfMotifs(dnaStrings)
# print(score)


# GREEDYMOTIFSEARCH(Dna, k,t)
# form a set of k-mers BestMotifs by selecting 1st k-mers in each string from Dna
# for each k-mer Motif in the 1st string from Dna
# Motif1 <- Motif
# for i = 2 to t
# form Profile from motifs Motif1, ., Motifi - 1
# Motifi <0 Profile-most probable k-mer in the i-th string in Dna
# Motifs <- (Motif1, ., Motift)
# if Score(Motifs) < Score(BestMotifs)
# BestMotifs <- Motifs
# output BestMotifs

GreedyMotifSearch <- function(dnaStrings, k, t, pseudocount=F){
  bestMotifs <- sapply(dnaStrings, substring, 1, k)
  bestMotifsScore <- scoreOfMotifs(bestMotifs)
  
  kmersInString1 <- getKmers(text=dnaStrings[1], k=k, unique=T)
  
  updateBestMotifsFromAnIteration <- function(initialMotif){
    curMotifs <- rep("", t)
    curMotifs[1] <- initialMotif
    
    tmp <- sapply(2:t, function(i){
      profile <- motifsToProfile(motifs=curMotifs[1:(i-1)], pseudocount=pseudocount)
      curMotifs[i] <<- profileMostProbableKmer(dna=dnaStrings[i], k=k, profile=profile)[1]
      return(NA)
    })
    
    curMotifsScore <- scoreOfMotifs(curMotifs)
    if(curMotifsScore < bestMotifsScore){
      bestMotifs <<- curMotifs
      bestMotifsScore <<- curMotifsScore
    }
    
  }
  #debug(updateBestMotifsFromAnIteration)
  
  sapply(kmersInString1, updateBestMotifsFromAnIteration)
  
  return(bestMotifs)
}

#debug(GreedyMotifSearch)

motifs <- GreedyMotifSearch(dnaStrings, k, t, pseudocount=T)
cat(motifs, sep="\n")