source("code/utility.motif.R")

inputFile <- "data/medium_string_data-3.txt"
inputs <- readLines(con=inputFile, warn=F)
numbers <- strsplit(inputs[1], " ")[[1]]
k <- as.integer(numbers[1])
dnaCollection <- inputs[2:length(inputs)]


# MEDIANSTRING(Dna, k)
# BestPattern <- AAA.AA
# for each k-mer Pattern from AAA.AA to TTT.TT
# if d(Pattern, Dna) < d(BestPattern, Dna)
# BestPattern <- Pattern
# output BestPattern

medianStringBruteForce <- function(dnaCollection, k){
  allPatterns <- getAllPossibleDna(k)
  distances <- sapply(allPatterns, distanceBetweenPatternAndDnaCollection, dnaCollection)
  minDistanceIndex <- which(distances==min(distances))
  medianString <- allPatterns[minDistanceIndex[1]]
  return(medianString)
}

medianStr <- medianStringBruteForce(dnaCollection, k)
cat(medianStr)