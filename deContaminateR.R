rm(list = ls())

library(stringr)
library(Biostrings)
library(dplyr)
library(readr)

#Code to fix yale sequences based on the contamination reports
#Set working directory
setwd("D://Documents/bmc_seq_submission")

#Function to read fasta file and make a data frame with sequence name and sequences
fasta4epitope <- function(fasta_file){
  #Read fasta file containing epitope sequences
  epitope_seq <- readAAStringSet(fasta_file) 
  #Save sequence names in a vector
  seq_name <- names(epitope_seq)
  #Save sequences in a vector
  sequences <- paste(epitope_seq)
  #Save these in a data frame
  seq_df <- data.frame(seq_name, sequences)
  colnames(seq_df) <- c("annotation","SEQUENCE")
  #Add a column with epitope name
  seq_df$source <- strsplit(fasta_file, "[.]")[[1]][1] 
  #Change all columns to character
  seq_df %>% mutate_all(as.character)
  #Add column for length of sequence
  seq_df$seq_length <- nchar(seq_df$SEQUENCE)
  #Return data frame
  return(seq_df)
}

#Make list of fasta files
file_list <- list.files(path = "fasta_files")
sample_list <- str_split_fixed(file_list, "contigs", 2)[,1]

#Make list of contamination reports
exclude_list <- list.files(path = "exclude")
trim_list <- list.files(path = "trim")

for(i in 1:length(sample_list)){
  #i <- 32
  print(paste("Starting analysis of ", sample_list[i], sep = ""))
  
  #Load fasta file
  fasta_filepath <- paste("fasta_files/",file_list[which(grepl(sample_list[i],file_list))],sep = "")
  fasta_df <- fasta4epitope(fasta_filepath)
  print(paste("Number of sequences WITH LENGTH LESS THAN 200: ", length(which(fasta_df$seq_length < 200)), sep = ""))
  print(paste("Number of sequences before exclusion: ", nrow(fasta_df), sep = ""))
  
  #Load exclude file
  exclude_filepath <- paste("exclude/",exclude_list[which(grepl(sample_list[i], exclude_list))],sep = "")
  exclude_file <- read_delim(exclude_filepath, delim = "\t")
  if(nrow(exclude_file != 0)){
    colnames(exclude_file)[1] <- "column_name" 
    exclude_nodelist <- str_split_fixed(exclude_file$column_name,"length_",2)[,1]
    #Exclude sequences from exclude nodes
    fasta_df <- subset(fasta_df, !grepl(paste(exclude_nodelist, collapse="|"), annotation))
  }
  print(paste("Number of sequences after exclusion: ", nrow(fasta_df), sep = ""))
  
  #Load trim file
  trim_filepath <- paste("trim/",trim_list[which(grepl(sample_list[i], trim_list))],sep = "")
  trim_file <- read_delim(trim_filepath, delim = "\t")
  if (nrow(trim_file) != 0){
    colnames(trim_file)[1] <- "column_name" 
    trim_file$node <- str_split_fixed(trim_file$column_name,"length_",2)[,1]
    trim_file$start <- str_split_fixed(trim_file$column_name,"\t",4)[,3]
    trim_file$end <- str_split_fixed(trim_file$start,"[.][.]",2)[,2]
    trim_file$start <- str_split_fixed(trim_file$start,"[.][.]",2)[,1]
    trim_file <- trim_file[,-1]
    
    for(j in 1:nrow(trim_file)){
      trim_seq <- substr(fasta_df$SEQUENCE[which(grepl(trim_file$node[j], fasta_df$annotation))], trim_file$start[j], trim_file$end[j])
      tmp_vec <- fasta_df$SEQUENCE[which(grepl(trim_file$node[j], fasta_df$annotation))]
      print(paste("Length of sequence at",trim_file$node[j], nchar(tmp_vec), sep = " "))
      fasta_df$SEQUENCE[which(grepl(trim_file$node[j], fasta_df$annotation))] <- gsub(trim_seq,"",tmp_vec)
      print(paste("Length of sequence after gsub", nchar(fasta_df$SEQUENCE[which(grepl(trim_file$node[j], fasta_df$annotation))]), sep = " "))
    }
  }
  
  fasta_df$seq_length <- nchar(fasta_df$SEQUENCE)
  #Drop sequences with length less than 200
  fasta_df <- subset(fasta_df, seq_length > 199)
  
  #Save results in output file
  output_filename <- paste("fasta_fixed/",file_list[i],sep = "")
  output_filename <- gsub("fsa","fasta",output_filename)
  write(paste(">",fasta_df$annotation,"\n",fasta_df$SEQUENCE,sep = ""),file = output_filename)
  
  print(paste("Completed analysis for ",  sample_list[i], sep = ""))
}

