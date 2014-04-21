#inputFile <- "data/6.2-manhattan-2.txt"
inputFile <- "C:/Users/Ashis/Downloads/dataset_73_9.txt"
outputFile <- "results/6.2-manhattan-1.txt"

inputs <- readLines(con=inputFile, warn=F)
numbers <- as.integer(strsplit(inputs[1], "\\s")[[1]])
row <- numbers[1] + 1
col <- numbers[2] + 1

downInput <- read.table(file=inputFile, header=F, quote="", colClasses="numeric",skip=1, nrows=row-1)
rightInput <- read.table(file=inputFile, header=F, quote="", colClasses="numeric",skip=1+row, nrows=row)

down <- matrix(data=0, nrow=row, ncol=col)
down[2:row,] <- as.matrix(downInput)
right <- matrix(data=0, nrow=row, ncol=col)
right[,2:col] <- as.matrix(rightInput)

len <- matrix(data=0, nrow=row, ncol=col)
len[,1] <- cumsum(down[,1])
len[1,] <- cumsum(right[1,])

for(j in 2:col){
  for(i in 2:row){
    len[i,j] <- max(len[i,j-1]+right[i,j], len[i-1,j]+down[i,j])
  }
}

 writeLines(paste(len[row,col]))

# DPChange <- function(money, coins){
#   minNumCoins <- 0
#   maxCoin <- max(coins)
#   for(m in 1:money){
#     candidates <- (m - coins) + 1    # +1 for indexing
#     if(m < maxCoin){
#       candidates <- candidates[candidates>0]
#     }
#     
#     minNumCoins[m+1] <- ifelse(length(candidates)>0, min(minNumCoins[candidates] + 1), Inf)
#   }
#   
#   return(minNumCoins[m+1])
# }
# 
# minCoins <- DPChange(money, coins)
# writeLines(paste(minCoins))
