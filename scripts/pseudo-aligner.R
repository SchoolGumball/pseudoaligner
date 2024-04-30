# this program is a basic breakdown of the methods used by kallisto in an attempt to better undestand the workings of the pseudo alligner
# https://arxiv.org/pdf/1505.02710.pdf
# plase note that the methods used for kallisto are more advanced, we just follow the big strokes with simpler methods 

#to do list
#read reference genome - DONE
#parse reference genome to a quick access data type - DONE
#load reads from file - DONE
#process reads to a quick access format - DONE
#use heuristics for fast alignment
#output the alignment somehow
#optional: test acurracy compared to kallisto (almost no point but might as well)

source("scripts/hash_sequence.R")
source("scripts/hash_reads.R")
source("scripts/aligner.R")

data = parse_fasta("raw_data/zika.fna")

index = create_hash_index(data)

#save the fasta so we dont need to reprocess later
saveRDS(data, file = "processed_data/data.rds")

#load the fasta from file (only relevant if we dont wanna run the index again)
index = readRDS("processed_data/data.rds")

#save the index so we dont need to reprocess later
saveRDS(index, file = "processed_data/index.rds")

#load the index from file (only relevant if we dont wanna run the index again)
index = readRDS("processed_data/index.rds")

#test read data set
#https://www.ncbi.nlm.nih.gov/sra/SRX7443658[accn] SRR10769568
#
#smaller data set for dev
#https://www.ncbi.nlm.nih.gov/sra/SRX7443590[accn] SRR10769636

reads = parse_fastq("raw_data/SRR10769636.fastq")

reads_index <- create_reads_index(reads)

#save the reads so we dont need to reprocess later
saveRDS(reads, file = "processed_data/reads.rds")

#load the reads from file (only relevant if we dont wanna run the reads index again)
reads = readRDS("processed_data/reads.rds")

#save the reads index so we dont need to reprocess later
saveRDS(reads_index, file = "processed_data/reads_index.rds")

#load the reads index from file (only relevant if we dont wanna run the reads index again)
reads_index = readRDS("processed_data/reads_index.rds")


#run the alignment
alignment <- align(index, reads_index)
