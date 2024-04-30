#hashes de reads

#in case it needs to access the class for the reads
source("scripts/read_fastq.R")

#facilitates oop because gumball is stubborn
library(R6)

#used to represent the reads index
readsIndex <- R6Class(
  "readsIndex",
  public = list(
    n_reads = NULL,
    accession_number = NULL,
    read_index = NULL,
    initialize = function(accession_number, n_reads, read_index) {
      self$accession_number <- accession_number
      self$n_reads <- n_reads
      self$read_index <- read_index
    },
    print = function() {
      cat("Accession Number(s):", paste(self$accession_number, collapse = " "), "\n")
      cat("Number of reads:", self$n_reads, "\n")
    }
  )
)

#trakes a readsData object returns an readsIndex object with the corresponding hash table, please don't break it
#the index will use k length mers, please dont use a crazy number, it will kill your cpu
#k MUST be the same as the reference index or this will die, might integrate it all later, no promises
create_reads_index <- function(reads_data, k = 10){
  #we use a list cause its efficient enough
  index = list()
  accession_number = reads_data$accession_number
  n_reads = reads_data$n_reads
  
  #access the reads
  r <- reads_data$read_data
  
  # Iterate through the reads
  for (i in names(r)){
    #we extract the fragment
    fragment <- r[[i]]
    
    kmer <- substr(fragment, start = 1, stop = k)
    
    #hash using  SHA-1, don't think too hard about this one
    hash_value = digest(kmer)
    
    # Check if the k-mer already exists in the hash table
    if (is.null(index[[i]])) {
      # If not, create a new list to store positions
      index[[i]] <- list()
    }
    
    # Add fragment id and sequence under the hash
    index[[i]] <- c(hash_value, fragment)
  }
  
  read_index = readsIndex$new(accession_number, n_reads, index)
  return (read_index)
}