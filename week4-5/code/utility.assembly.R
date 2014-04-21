library(graph)
library(data.table)

source("code/utility.motif.R")
source("code/utility.graph.R")

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


getDeBrujinGraph <- function(patterns){
  #patterns <- sort(unique(patterns))
  k <- str_length(patterns[1])
  
#   from <- substr(patterns, 1, k-1)
#   to <- substr(patterns, 2, k)
#   allEdges <- data.frame(from=from, to=to)
#   allEdges <- data.table(allEdges)

  allEdges <- data.table(from=substr(patterns, 1, k-1), to=substr(patterns, 2, k))
  
  #collapsedEdges <- aggregate(to~from, allEdges, paste0, sep="", collapse=",")
  collapsedEdges <- allEdges[, paste0(to, sep="", collapse=","), by = from]
  #collapsedEdges <- as.data.frame(collapsedEdges, stringsAsFactors=F)
  setnames(collapsedEdges,old=c("from", "to"))
  #colnames(collapsedEdges) <- c("from","to")
  return(collapsedEdges)
  
#   patterns <- sort(unique(patterns))
#   names(patterns) <- patterns
#   
#   #edgesFrom <- c()
#   #edgesTo <- c()
#   
#   createDeBrujinEdge <- function(pattern){
#     len <- str_length(pattern)
#     prefix <- substr(pattern,start=1,stop=len-1)
#     suffix <- substr(pattern,start=2,stop=len)
#     return(c(from=prefix, to=suffix))
#   }
#   
#   allEdges <- sapply(patterns, createDeBrujinEdge)
#   allEdges <- t(allEdges)
#   collapsedEdges <- aggregate(to~from, allEdges, paste0, sep="", collapse=",")
#   return(collapsedEdges)

}

getDeBrujinGraphFromPairedReads <- function(patterns, sep=NA){
  k <- (str_length(patterns[1])-1)/2
  if(is.na(sep) || is.null(sep)){
    sep <- substr(patterns[1],k+1,k+1)
  }
  
  from <- paste(substr(patterns, 1, k-1), substr(patterns, k+2, 2*k), sep=sep)
  to <- paste(substr(patterns, 2, k), substr(patterns, k+3, 2*k+1), sep=sep)
  allEdges <- data.table(from=from, to=to)
  collapsedEdges <- allEdges[, paste0(to, sep="", collapse=","), by = from]
  setnames(collapsedEdges,old=c("from", "to"))
  return(collapsedEdges)
}

getDeBrujinGraph_v1 <- function(patterns){
  
  patterns <- sort(unique(patterns))
  names(patterns) <- patterns
  
  #edgesFrom <- c()
  #edgesTo <- c()
  
  createDeBrujinEdge <- function(pattern){
    len <- str_length(pattern)
    prefix <- substr(pattern,start=1,stop=len-1)
    suffix <- substr(pattern,start=2,stop=len)
    return(c(from=prefix, to=suffix))
  }
  
  allEdges <- sapply(patterns, createDeBrujinEdge)
  allEdges <- t(allEdges)
  collapsedEdges <- aggregate(to~from, allEdges, paste0, sep="", collapse=",")
  return(collapsedEdges)
}


writeEdgesInFile <- function(from, to, file, sep=" -> "){
  write.table(x=data.frame(from=from, to=to), file=file, quote=F, sep=sep, row.names=F, col.names=F)
  return(NA)
}

## each node is a k-mer.
genomeFromEulerPath <- function(epath){
  k <- str_length(epath[1])
  lastBits <- sapply(epath,substr,k,k)
  genome <- paste0(substr(epath[1],1,k-1), paste0(lastBits,collapse=""))
  return(genome)
}

## each node is a paired k-mer.
genomeFromEulerPathOfPairReads <- function(epath, d){
  k <- (str_length(epath[1])-1)/2
  startBits1 <- substr(epath[1],1,k-1)
  lastBits1 <-paste0(substr(epath,k,k), collapse="")
  lastBits2 <- paste0(substr(epath[(length(epath)-k-d):length(epath)],2*k+1,2*k+1), collapse="")
  genome <- paste0(startBits1, lastBits1, lastBits2)
  return(genome)
}

# g = debrujin graph
# g has to be connected
contigs <- function(g){
  nodes <- nodes(g)
  degs <- degree(g)
  
  startFlag <- (degs$outDegree >=1) &  !(degs$inDegree==1 & degs$outDegree==1)
  stopFlag <- degs$outDegree!=1 | degs$inDegree>1 
  
  #handle cycle with 1 outgoing and 1 incoming edge
  if(all(!startFlag)){
    startFlag[1] = TRUE
    stopFlag[1] = TRUE
  }
  
  ed <- edges(g)
  
  startNodes <- rep(nodes[startFlag], degs$outDegree[startFlag])
  contigs <- character()
  
  for(i in 1:length(startNodes)){
    walk <- startNodes[i]
    curNode <- startNodes[i]
    
    while(TRUE){
      nextNode <- ed[[curNode]][1]
      walk <- append(walk, nextNode)
      ed[[curNode]] <- ed[[curNode]][-1]
      
      if(stopFlag[nextNode]){
        break;
      }
      
      curNode <- nextNode
    }
    
    contigs[i] <- genomeFromEulerPath(epath=walk)
  }
  
  
  return(contigs)
  
}



