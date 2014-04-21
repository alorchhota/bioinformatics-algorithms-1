# Approximate Pattern Matching Problem: Find all approximate occurrences of a pattern in a string.
# Input: Two strings Pattern and Text along with an integer d.
# Output: All positions where Pattern appears in Text with at most d mismatches.

library("stringr")

patternMatchApproximate <- function(string, pattern, d, indexBase=1){
  
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

# pattern <- "ATTCTGGA"
# text <- "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAAT"
# d <- 3
# patternMatchApproximate(text, pattern, d, indexBase=0)


inputFile <- "data/approximate_match_data_3.txt"
inputs <- readLines(con=inputFile, n=3)
pattern <- inputs[1]
text <- inputs[2]
d <- as.integer(str_trim(inputs[3]))
matchIndexes <- patternMatchApproximate(text, pattern, d, indexBase=0)
cat(matchIndexes)

