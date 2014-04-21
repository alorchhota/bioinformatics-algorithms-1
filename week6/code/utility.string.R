require(stringr)
lcs <- function(s,t){
  s <- strsplit(s,split="")[[1]]
  t <- strsplit(t,split="")[[1]]
  
  row <- length(s) + 1
  col <- length(t) + 1
  
  len <- matrix(data=0, nrow=row, ncol=col)
  len[,1] <- 0
  len[1,] <- 0
  
  # Backtrack: 1:Down, 2:Right, 3:Diagonal
  backtrack <- matrix(data=0, nrow=row, ncol=col)
  backtrack[,1] <- 1
  backtrack[1,] <- 2
  
#   rownames(backtrack) <- c("-",s)
#   colnames(backtrack) <- c("-",t)
#   rownames(len) <- c("-",s)
#   colnames(len) <- c("-",t)
  
  for(j in 2:col){
    for(i in 2:row){
      candidates <- c(len[i-1,j], len[i,j-1], len[i-1,j-1]+as.integer(s[i-1]==t[j-1]))
      m <- which.max(candidates)
      len[i,j] <- candidates[m]
      backtrack[i,j] <- m
    }
  }
  
  lcslen <- len[row,col]
  lcs <- character(length=lcslen)
  i <- row
  j <- col
  
  while(lcslen>0){
    if(backtrack[i,j]==3){
      lcs[lcslen] <- s[i-1]
      lcslen <- lcslen-1
      i <- i-1
      j <- j-1
    } else if(backtrack[i,j]==2){
      j <- j-1
    } else if(backtrack[i,j]==1){
      i <- i-1
    }
  }
  
  lcs <- paste0(lcs, collapse="")
  return(lcs)
}