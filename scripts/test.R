#this is a test for the aligner jsut run this, trust me

source("scripts/hash_sequence.R")
source("scripts/hash_reads.R")
source("scripts/aligner.R")

#light test, runs really quick 

#process index
data = parse_fasta("raw_data/zika.fna")

index = create_hash_index(data, 5)

#process reads
reads = parse_fastq("raw_data/SRR10769636.fastq")

reads_index <- create_reads_index(reads, 5)

#run the alignment
alignment <- align(index, reads_index)

#this is the actual data from the alignment as in the sequences aligned and where they go
position_data <- alignment$alignment_data

######
######WAY HEAVIER TEST PLEASE DO RUN THIS BUT KEEP IN MIND IT CAN POTENTIALLY RUN FOR A WHILE
######LIKE FOR REAL THIS ISN'T EVEN THAT BIG AN ALIGNMENT FILE AND IT'S STILL HUGE!!!!!
######

#process index
data = parse_fasta("raw_data/zika.fna")

index = create_hash_index(data, 5)

#process reads
reads = parse_fastq("raw_data/SRR10769568.fastq")

reads_index <- create_reads_index(reads, 5)

#run the alignment
alignment <- align(index, reads_index)

#this is the actual data from the alignment as in the sequences aligned and where they go
position_data <- alignment$alignment_data