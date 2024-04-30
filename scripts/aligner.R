#this takes care of the alignments, best of luck

#facilitates oop because gumball is stubborn
library(R6)
#for the df 
library(dplyr)

#used to represent the reads index
Alignment <- R6Class(
  "Alignment",
  public = list(
    reference_accession_number = NULL,
    reference_common_name = NULL,
    n_reads = NULL,
    reads_accession_number = NULL,
    alignment_data = NULL,
    unalligned_data = NULL,
    reference_sequence = NULL,
    initialize = function(reference_accession_number,  reference_common_name, n_reads, reads_accession_number, alignment_data, unalligned_data, reference_sequence) {
      self$reference_accession_number <- reference_accession_number
      self$reference_common_name <- reference_common_name
      self$n_reads <- n_reads
      self$reads_accession_number <- reads_accession_number
      self$alignment_data <- alignment_data
      self$unalligned_data <- unalligned_data
      self$reference_sequence <- reference_sequence
    },
    print = function() {
      cat("Reference accession Number:", self$reference_accession_number, "\n")
      cat("Common Name:", self$reference_common_name, "\n")
      cat("Number of reads:", self$n_reads, "\n")
      cat("Reads accession Number(s):", paste(self$reads_accession_number, collapse = " "), "\n")
    }
  )
)

#this will do the alignment effort, fingers crossed it isn't supper heavy
align <- function(index, reads_index) {
  r_index <- reads_index$read_index
  i <- index$hash_index
  alignment <-  data.frame(fragment_id = character(),
                           fragment_sequence = character(),
                           aligned_positions = list())
  unaligned <- list()
  
  for (n in names(r_index)){
    #we get the hash
    h <- r_index[[n]][1]
    #we check if the index matches any entries
    if (!is.null(i[[h]])){
      #we extract the info we want
      #the order is dependant on the readsIndex data type in case this breaks, probably why
      id <- n
      frag <- r_index[[n]][2]
      pos <- i[[h]]
      
      #then we add to df
      #https://stackoverflow.com/questions/28467068/how-to-add-a-row-to-a-data-frame-in-r
      alignment <- bind_rows(alignment, data.frame(fragment_id = id,
                                                   fragment_sequence = frag,
                                                   aligned_positions = unlist(pos)))
      
    } else{ #if not we throw it into the unaligned list
      unaligned[[n]] <- r_index[[n]][2]
    }
    
  }
  final_alignment <- Alignment$new(index$accession_number, index$common_name, reads_index$n_reads, reads_index$accession_number, alignment, unaligned, index$seq)
  return(final_alignment)
}