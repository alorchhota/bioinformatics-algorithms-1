inputFile <- "data/change_problem-2.txt"
inputFile <- "C:/Users/Ashis/Downloads/dataset_72_9.txt"
outputFile <- "results/6.1-change-problem.txt"

inputs <- readLines(con=inputFile, warn=F)
money <- as.integer(inputs[1])
numbers <- strsplit(inputs[2], ",")[[1]]
coins <- as.integer(numbers)


DPChange <- function(money, coins){
  minNumCoins <- 0
  maxCoin <- max(coins)
  for(m in 1:money){
    candidates <- (m - coins) + 1    # +1 for indexing
    if(m < maxCoin){
      candidates <- candidates[candidates>0]
    }
    
    minNumCoins[m+1] <- ifelse(length(candidates)>0, min(minNumCoins[candidates] + 1), Inf)
  }
  
  return(minNumCoins[m+1])
}

minCoins <- DPChange(money, coins)
writeLines(paste(minCoins))
