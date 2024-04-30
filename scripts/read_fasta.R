#takes care of reading and parsing a fasta file into a class for us to work

#facilitates oop because gumball is stubborn
#https://r6.r-lib.org/articles/Introduction.html
library(R6)

# Define a class to hold genome data
GenomeData <- R6Class(
  "GenomeData",
  public = list(
    accession_number = NULL,
    common_name = NULL,
    seq = NULL,
    initialize = function(accession_number, common_name, seq) {
      self$accession_number <- accession_number
      self$common_name <- common_name
      self$seq <- seq
    },
    print = function() {
      cat("Accession Number:", self$accession_number, "\n")
      cat("Common Name:", self$common_name, "\n")
    }
  )
)

# Function to parse FASTA file and return a GenomeData object
# will break if FASTA has multiple sequences, will only give back last one, don't ask me to change it
parse_fasta <- function(file_path) {
  genome_data = NULL
  accession_numer = NULL
  name = NULL
  seq = NULL
  
  # Read FASTA file line by line
  lines <- readLines(file_path)
  for (line in lines) {
    #is its the header
    if (startsWith(line, ">")) {
      # Extract accession number and common name from header line
      # we first find where the information is in the line
      #https://stackoverflow.com/questions/14249562/find-the-location-of-a-character-in-string
      name_start_index = unlist(lapply(gregexpr(pattern = ' ', line), min)) + 1
      name_end_index = unlist(lapply(gregexpr(pattern = ',', line), min)) - 1
      
      #substring using the found index
      accession_numer = substr(line, 2, name_start_index - 2)
      name = substr(line, name_start_index, name_end_index)
    }
    #anything else, hope its an actual sequence
    else{
      seq = paste(seq, line, sep = "")
    }
    
  }
  #structure and return out data
  genome_data = GenomeData$new(accession_numer, name, seq)
  return(genome_data)
}