library(igraph)
library(graph)


# Graph File Format:
# 0 -> 3
# 1 -> 0
# 2 -> 1,6
# 3 -> 2
# 4 -> 2
# 5 -> 4
# 6 -> 5,8
# 7 -> 9
# 8 -> 7
# 9 -> 6
  
graphFromFile_v0 <- function(inputFile){
  inputs <- readLines(con=inputFile, warn=F)
  splittedInputs <- strsplit(inputs, "\\s->\\s")
  elems <- unlist(splittedInputs)
  len <- length(elems)
  from <- elems[seq(from=1, to=len, by=2)]
  to <- elems[seq(from=2, to=len, by=2)]
  edL <- strsplit(to, ",")
  #edL <- lapply(edL, function(l){return(as.numeric(l))}) ## convert to numeric edge list
  names(edL) <- from
  gR <- graphNEL(nodes=from, edgeL=edL, edgemode="directed")
  return(gR)
}

graphFromFile <- function(inputFile){
  inputs <- readLines(con=inputFile, warn=F)
  splittedInputs <- strsplit(inputs, "\\s->\\s")
  elems <- unlist(splittedInputs)
  len <- length(elems)
  from <- elems[seq(from=1, to=len, by=2)]
  to <- elems[seq(from=2, to=len, by=2)]
  edL <- strsplit(to, ",")
  names(edL) <- from
  
  ## add non existing nodes and correspoding empty edge
  noEdgeNodes <- setdiff(unique(unlist(edL)), from)
  from <- c(from, noEdgeNodes)
  lapply(noEdgeNodes, function(n){edL[[n]] <<- new("character")})
  
  gR <- graphNEL(nodes=from, edgeL=edL, edgemode="directed")
  return(gR)
}

# Format: 
# 0->1:7
# 0->2:4
# 2->3:2
# 1->4:1
# 3->4:3  

#each edge is a string

weightedGraphFromEdges <- function(edges, sep.dir="->", sep.w=":"){
  splittedInputs <- strsplit(edges, split=paste0(sep.dir,"|",sep.w))
  elems <- unlist(splittedInputs)
  len <- length(elems)
  from <- elems[seq(from=1, to=len, by=3)]
  to <- elems[seq(from=2, to=len, by=3)]
  w <- elems[seq(from=3, to=len, by=3)]
  
  nodes <- unique(append(from, to))
  g <- graphNEL(nodes=nodes, edgemode="directed")
  g <- addEdge(from=from, to=to, weights=w, graph=g)
  return(g)
}

eulerCycle <- function(g, start=NULL){
  eulerCycle <- c()
  curNode <- ifelse(is.null(start), nodes(g)[1], start)
  
  
  while(!is.na(curNode)){
    cycle <- curNode
    
    while(!is.na(nextNode <- randomWalkNext(g, curNode))){
      g <- removeEdge(graph=g, from=curNode, to=nextNode)
      cycle <- append(cycle, nextNode)
      curNode <- nextNode
    }
    
    if(length(eulerCycle)==0){
      eulerCycle <- cycle
    }
    else{
      insertIndex <- which(eulerCycle==cycle[1])[1]
      eulerCycle <- append(eulerCycle,after=insertIndex,values=cycle[-1])
    }
    
    curNode <- getAnUnexploredNode(g, nodes=eulerCycle)
  }
  
  return(eulerCycle)
  
}

getAnUnexploredNode <- function(g, nodes){
  degrees <- degree(g, Nodes=nodes)
  nodeIndexes <- which(degrees$outDegree+degrees$inDegree>0)
  node <- NA
  if(length(nodeIndexes)>0){
    node <- nodes[nodeIndexes[1]]
  }
  return(node)
}

getUnexploredNodes <- function(g, nodes){
  degrees <- degree(g, Nodes=nodes)
  nodeIndexes <- which(degrees$outDegree>0)
  return(nodes[nodeIndexes])
}

randomWalkNext <- function(g, from){
  outEdges <- edges(object=g, which=from)[[1]]
  nextNode <- NA
  if(length(outEdges)>0){
    nextNode <- outEdges[1]
  }
  return(nextNode)
}


eulerian <- function(g, start=NULL){
  eulerPath <- c()
  
  curNode <- ifelse(is.null(start), selectAStartNodeForEulerian(g), start)
  
  unexplored <- character()
  count <- 0
  
  ed <- edges(g)
  startIndex <- 1
  endIndex <- length(ed)
  
  while(!is.na(curNode)){
    cycle <- curNode
    
    while(TRUE){
      
      # !is.na(nextNode <- randomWalkNext(g, curNode))
      moves <- ed[[curNode]]
      nextNode <- ifelse(length(moves)>0, moves[1], NA)
      if(is.na(nextNode))
        break;
      
      #g <- removeEdge(graph=g, from=curNode, to=nextNode)
      ed[[curNode]] <- ed[[curNode]][-1]
      
      
      cycle <- append(cycle, nextNode)
      curNode <- nextNode
      
      count <- count + 1
      if(count%%1000==0){
        print(count)
      }
      
      
    }
    
    if(length(eulerPath)==0){
      eulerPath <- cycle
    }
    else{
      insertIndex <- which(eulerPath==cycle[1])[1]
      eulerPath <- append(eulerPath,after=insertIndex,values=cycle[-1])
    }
    
    
    #t1 <- Sys.time()
    curNode <- NA
    while((startIndex <- startIndex + 1) <= endIndex){
      node <- eulerPath[startIndex]
      # if has unexplored edge
      if(length(ed[[node]])>0){
        curNode <- node
        break
      }
    }
    
    #print(curNode)
    
    #t2 <- Sys.time()
    #print(t2 - t1)
    
#     t1 <- Sys.time()
#     
#     hasEdge <- sapply(eulerPath, function(n){
#       return(length(ed[[n]])>0)
#     })
#     
#     indexes <- which(hasEdge)
#     curNode <- ifelse(length(indexes)>0, eulerPath[indexes[1]], NA)
#     
#     t2 <- Sys.time()
#     print(t2 - t1)
    
    count <- count + 1
    #     if(count%%1000==0){
    #       print(count)
    #     }
  }
  
  #print(count)
  return(eulerPath)
  
}


eulerian_v2 <- function(g, start=NULL){
  eulerPath <- c()
  
  curNode <- ifelse(is.null(start), selectAStartNodeForEulerian(g), start)
  count <- 0
  
  ed <- edges(g)
  
  while(!is.na(curNode)){
    cycle <- curNode
    
    while(TRUE){
      
      # !is.na(nextNode <- randomWalkNext(g, curNode))
      moves <- ed[[curNode]]
      nextNode <- ifelse(length(moves)>0, moves[1], NA)
      if(is.na(nextNode))
        break;
      
      #g <- removeEdge(graph=g, from=curNode, to=nextNode)
      ed[[curNode]] <- ed[[curNode]][-1]
      
      
      cycle <- append(cycle, nextNode)
      curNode <- nextNode
      
#       count <- count + 1
#       if(count%%1000==0){
#         print(count)
#         #gc()
#       }
        
        
    }
    
    if(length(eulerPath)==0){
      eulerPath <- cycle
    }
    else{
      insertIndex <- which(eulerPath==cycle[1])[1]
      eulerPath <- append(eulerPath,after=insertIndex,values=cycle[-1])
    }
    
    #unexplored <- getUnexploredNodes(g, nodes=unique(append(unexplored, cycle)))
    #curNode <- ifelse(length(unexplored)==0, NA, unexplored[1])
    
    #g <- graphNEL(nodes=nodes(g), edgeL=ed, edgemode="directed")
    #curNode <- getAnUnexploredNode(g, nodes=eulerPath)

    t1 <- Sys.time()

    hasEdge <- sapply(eulerPath, function(n){
      return(length(ed[[n]])>0)
    })

    indexes <- which(hasEdge)
    curNode <- ifelse(length(indexes)>0, eulerPath[indexes[1]], NA)

    t2 <- Sys.time()
    print(t2 - t1)

    count <- count + 1
#     if(count%%1000==0){
#       print(count)
#     }
  }
  
  print(count)
  return(eulerPath)
  
}


eulerian_v1 <- function(g, start=NULL){
  eulerPath <- c()
  
  curNode <- ifelse(is.null(start), selectAStartNodeForEulerian(g), start)
  
  
  while(!is.na(curNode)){
    cycle <- curNode
    
    while(!is.na(nextNode <- randomWalkNext(g, curNode))){
      g <- removeEdge(graph=g, from=curNode, to=nextNode)
      cycle <- append(cycle, nextNode)
      curNode <- nextNode
    }
    
    if(length(eulerPath)==0){
      eulerPath <- cycle
    }
    else{
      insertIndex <- which(eulerPath==cycle[1])[1]
      eulerPath <- append(eulerPath,after=insertIndex,values=cycle[-1])
    }
    
    curNode <- getAnUnexploredNode(g, nodes=eulerPath)
  }
  
  return(eulerPath)
  
}


selectAStartNodeForEulerian <- function(g){
  nodes <- nodes(g)
  deg <- degree(g, Nodes=nodes)
  degDiff <- deg$outDegree - deg$inDegree
  degDiff1Index <- which(degDiff==1)
  if(length(degDiff1Index)==1){
    startNode <- nodes[degDiff1Index]
  }else if(length(degDiff1Index)==0){
    startNode <- nodes[1]
  }else{
    stop("There is no eulerian path in this graph.")
  }
  
  return(startNode)
}


# getDeBrujinEdges <- function(patterns){
#   k <- length(patterns[1])
#   from <- substr(patterns, 1, k-1)
#   to <- substr(patterns, 2, k)
#     
#   substr(c("abc","def"),1,2)
#   
#   edges <- lapply(patterns, function(p){return(character())})
#   names(edges) <- patterns
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
# }

