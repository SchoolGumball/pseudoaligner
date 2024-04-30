#this one reads the alignments 

#facilitates oop because gumball is stubborn
library(R6)

readsData <- R6Class(
  "readsData",
  public = list(
    n_reads = NULL,
    accession_number = NULL,
    read_data = NULL,
    initialize = function(accession_number, read_data, n_reads) {
      self$accession_number <- accession_number
      self$read_data <- read_data
      self$n_reads <- n_reads
    },
    print = function() {
      cat("Accession Number(s):", paste(self$accession_number, collapse = " "), "\n")
      cat("Number of reads:", self$n_reads, "\n")
    }
  )
)

# Function to parse FASTQ file and return a readsData object
# will break if header has more "." than expected again don't @ me, please just don't break it
parse_fastq <- function(file_path) {
  reads_data = NULL
  accession_numer = list()
  rd = list()
  save_next = FALSE
  c_id = NULL
  n_reads = 0
  
  # Read FASTQ file line by line
  lines <- readLines(file_path) #wonder if this will break
  for (line in lines) {
    #is its the start of a sequence
    if (startsWith(line, "@")) {
      # Extract accession number and 
      # we first find where the information is in the line
      #https://stackoverflow.com/questions/14249562/find-the-location-of-a-character-in-string
      acc_end_index = unlist(lapply(gregexpr(pattern = '\\.', line), min)) - 1
      id_end_index = unlist(lapply(gregexpr(pattern = '\\.', line), max)) - 1
      
      #substring using the found index
      acc_numer = substr(line, 2, acc_end_index - 1)
      c_id = substr(line, 2, id_end_index)
      
      #find if acc number is new
      if (!(acc_numer %in% accession_numer)){
        #if not we add it
        accession_numer <- append(accession_numer, acc_numer)
      }
      
      #we then save the full id with a mock value
      rd[[c_id]] <- "temp"
      
      #we then note we have to save next line
      save_next = TRUE
      
      #jump a cicle
      next
    }
    
    #we save next line if relevant
    if (save_next){
      rd[[c_id]] <- line
      
      #do not forget this
      save_next = FALSE
      
      #jank wat to count reads
      n_reads <- n_reads + 1
    }
    
  }
  
  #structure and return out data
  reads_data = readsData$new(accession_numer, rd, n_reads)
  return(reads_data)
}