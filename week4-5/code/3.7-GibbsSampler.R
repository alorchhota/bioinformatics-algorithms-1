source("code/utility.motif.R")

#inputFile <- "data/gibbs_sampler-1.txt"
inputFile <- "C:/Users/Ashis/Downloads/dataset_43_4.txt"
inputs <- readLines(con=inputFile, warn=F)
numbers <- strsplit(inputs[1], " ")[[1]]
k <- as.integer(numbers[1])
t <- as.integer(numbers[2])
N <- as.integer(numbers[3])
dnaStrings <- inputs[2:length(inputs)]

#N <- min(N, 1200)

GibbsSampler <- function(dnaStrings, k, t, N, iteration=1){
  
  selectRandomKmer <- function(text, k){
    maxStart <- str_length(text) - k + 1
    randomKmerStart <- sample(1:maxStart, size=1)
    randomKmer <- substr(x=text, start=randomKmerStart, stop=randomKmerStart+k-1)
    return(randomKmer)
  }
  
  dnaStringsKmers <- lapply(dnaStrings, getKmers, k=k, unique=F)
  
  oneIterationOfGibbsSampler <- function(){
    curMotifs <- sapply(dnaStrings, selectRandomKmer, k)
    bestMotifs <- curMotifs
    bestMotifscountMatrix <- motifsToCountMatrix(bestMotifs)
    bestMotifsScore <- scoreOfMotifs(bestMotifs)
    curCountMatrix <- bestMotifscountMatrix
    
    coreLoopInGibbsSampler <- function(n){
      i <- sample(1:t, size=1)
      
      deletedMotif <- bestMotifs[i]
      selectedMotifs <- bestMotifs[-i]
      selectedCountMatrix <-  updateCountMatrixForMotifDeletion(countMatrix=bestMotifscountMatrix, motif=deletedMotif)
      selectedProfile <- countMatrixToProfile(selectedCountMatrix, pseudocount=T)
      curMotifs <- bestMotifs
      
      motif_i <- profileRandomKmer(profile=selectedProfile, kmers=dnaStringsKmers[[i]])
      curCountMatrix <- updateCountMatrixForMotifAddition(countMatrix=selectedCountMatrix, motif=motif_i)
      curMotifs[i] <- motif_i
      
      curMotifsScore <- scoreOfMotifsFromCountMatrix(countMatrix=curCountMatrix)
      if(curMotifsScore < bestMotifsScore){
        bestMotifs <<- curMotifs
        bestMotifsScore <<- curMotifsScore
        bestMotifscountMatrix <<- curCountMatrix
      }
      return(NA)
    }
    
    lapply(1:N, coreLoopInGibbsSampler)
    cat(file="out.txt", bestMotifs, bestMotifsScore, "===============", sep="\n", append=T)
    return(list(motifs=bestMotifs, score=bestMotifsScore))
    
  }
  
  vars <- ls(envir=globalenv())
  cl <- makeCluster(getOption("cl.cores", 4))
  clusterEvalQ(cl,{
    source("code/utility.motif.R")
  })
  
  clusterExport(cl=cl, varlist=vars)
  iteratedMotifs <-  parSapply(cl, 1:iteration, function(iter){
    oneIterationBestMotifs <- oneIterationOfGibbsSampler()
    return(oneIterationBestMotifs)
  })
  
  stopCluster(cl)
  
  scores <- unlist(iteratedMotifs["score",])
  bestIndex <- which(scores==min(scores))[1]
  bestMotifs <- iteratedMotifs["motifs",bestIndex]$motifs
  
  return(bestMotifs)
  
}


print(Sys.time())
motifs <- GibbsSampler(dnaStrings, k, t, N, iteration=20)
print(Sys.time())
cat(motifs, sep="\n")
