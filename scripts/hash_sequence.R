#takes care of creating the hash table

#contains the "class" and function to create the index table
source("scripts/read_fasta.R")

#facilitates oop because gumball is stubborn
library(digest)

#used to represent the hash index
Index <- R6Class(
  "Index",
  public = list(
    accession_number = NULL,
    common_name = NULL,
    seq = NULL,
    hash_index = NULL,
    k = NULL,
    initialize = function(accession_number, common_name, seq, hash_index, k) {
      self$accession_number <- accession_number
      self$common_name <- common_name
      self$seq <- seq
      self$hash_index <- hash_index
      self$k <- k
    },
    print = function() {
      cat("Accession Number:", self$accession_number, "\n")
      cat("Common Name:", self$common_name, "\n")
      cat("K-mers of length:", self$k, "\n")
    }
  )
)

#trakes a GenomeData object returns an Index object with the corresponding hash table, please don't break it
#the index will use k length mers, please dont use a crazy number, it will kill your cpu
create_hash_index <- function(genome_data, k = 10){
  #we use a list cause its efficient enough
  index = list()
  accession_number = genome_data$accession_number
  common_name = genome_data$common_name
  seq = genome_data$seq
  
  # Iterate through the sequence, rolling window
  for (i in 1:(nchar(seq) - k + 1)) {
    # Extract the k-mer
    kmer <- substr(seq, i, i + k - 1)
    
    #hash using  SHA-1, don't think too hard about this one
    hash_value = digest(kmer)
    
    # Check if the k-mer already exists in the hash table
    if (is.null(index[[hash_value]])) {
      # If not, create a new list to store positions
      index[[hash_value]] <- list()
    }
    
    # Add the position of the k-mer to the list
    # the c() is used to deal with a k-mer apearing in multiple places
    index[[hash_value]] <- c(index[[hash_value]], i)
  }
  
  full_index = Index$new(accession_number, common_name, seq, index, k)
  return (full_index)
}