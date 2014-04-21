# Frequent Words Problem: Find the most frequent k-mers in a string.
# Input: A string Text and an integer k.
# Output: All most frequent k-mers in Text.

library("stringr")

freuent.k.mers <- function(text, k){
  kmers <- sapply(1:(str_length(text)-k+1), function(start){
    end <- start + k -1
    si <- substr(x=text, start=start, stop=end)
    return(si)
  })
  
  freq <- table(kmers)
  maxFreq <- max(freq)
  maxFreqItems <- names(which(freq==maxFreq))
  return(list(frequent.k.mers=maxFreqItems, max.frequency=maxFreq))
}

#debug(freuent.k.mers)

text <- "ACGTTGCATGTCGCATGATGCATGAGAGCT"
k <- 4
kMers <- freuent.k.mers(text, k)
cat(kMers$frequent.k.mers)