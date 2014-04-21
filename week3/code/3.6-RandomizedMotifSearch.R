source("code/utility.motif.R")

inputFile <- "data/randomized-6.txt"
inputs <- readLines(con=inputFile, warn=F)
numbers <- strsplit(inputs[1], " ")[[1]]
k <- as.integer(numbers[1])
t <- as.integer(numbers[2])
dnaStrings <- inputs[2:length(inputs)]


# profile <- motifsToProfile(dnaStrings)
# score <- scoreOfMotifs(dnaStrings)
# print(score)


# RANDOMIZEDMOTIFSEARCH(Dna, k, t)
# randomly select k-mers Motifs = (Motif1, ., Motift) in each string from Dna
# BestMotifs <- Motifs
# while forever
# Profile <- Profile(Motifs)
# Motifs <- Motifs(Profile, Dna)
# if Score(Motifs) < Score(BestMotifs)
# BestMotifs <- Motifs
# else
#   output BestMotifs
# return

RandomizedMotifSearch <- function(dnaStrings, k, t, iteration=1){
  
  selectRandomKmer <- function(text, k){
    kmers <- getKmers(text=text, k=k, unique=T)
    randomKmer <- sample(x=kmers, size=1)
    return(randomKmer)
  }
  
  oneIterationOfRandomizedSearch <- function(){
    bestMotifs <- sapply(dnaStrings, selectRandomKmer, k)
    bestMotifsScore <- scoreOfMotifs(bestMotifs)
    
        
#     count <- 0
    while(TRUE){
#       count <- count + 1
#       cat (count)
#       cat(",")
      profile <- motifsToProfile(bestMotifs, pseudocount=T)
      newMotifs <- profileToMotifs(profile=profile, dnaStrings=dnaStrings)
      
      newMotifsScore <- scoreOfMotifs(newMotifs)
      if(newMotifsScore < bestMotifsScore){
        bestMotifs <- newMotifs
        bestMotifsScore <- newMotifsScore
      }
      else{
        return(list(motifs=bestMotifs, score=bestMotifsScore))
      }
      
    }
  }
  
  iteratedMotifs <- sapply(1:iteration, function(iter){
    if(iter %% 10 == 0)
      print(iter)
    oneIterationBestMotifs <- oneIterationOfRandomizedSearch()
    return(oneIterationBestMotifs)
  })
  
  scores <- unlist(iteratedMotifs["score",])
  bestIndex <- which(scores==min(scores))[1]
  bestMotifs <- iteratedMotifs["motifs",bestIndex]$motifs
  
  return(bestMotifs)
}


#debug(RandomizedMotifSearch)

motifs <- RandomizedMotifSearch(dnaStrings, k, t, iteration=100)
cat(motifs, sep="\n")
