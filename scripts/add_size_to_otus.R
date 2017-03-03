#!/usr/bin/env Rscript

# This script add to the sequence header its count from an OTU table.
# See https://github.com/alex-bagnoud/ArsM-evolutionary-placements/ for more details about this script.

library("optparse")

option_list = list(
    make_option(c("-t", "--otuTable"), type="character", default="1-reid_otus/11-prot_otu_table.txt", 
                help="OTU table input file", metavar="character"),
    make_option(c("-c", "--countColumn"), type="numeric", default=2, 
                help="Column number that contains count information in the OTU table file", metavar="numeric"),
    make_option(c("-i", "--fastaInput"), type="character", default="1-reid_otus/9-prot_otus.fasta", 
                help="Fasta input file", metavar="character"),
    make_option(c("-f", "--fastaOutput"), type="character", default="1-reid_otus/12-prot_otus_size.fasta", 
                help="Output fasta file name [default= %default]", metavar="character")
    ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# Import OTU table
otu_table_df <- read.table(opt$otuTable, header = FALSE)

# Import fasta file
library("Biostrings")
prot_string <- readAAStringSet(opt$fastaInput)
prot_df <- data.frame("headers" = names(prot_string), "seq" = paste(prot_string))

# Merge both data frames
merged_df <- merge(otu_table_df, prot_df, by.x = "V1", by.y = "headers")

# Create new headers
merged_df$new_headers <- paste(merged_df$V1, ";size=", merged_df[,opt$countColumn],";", sep = "")

# Export the new fasta file
new_fasta_string = AAStringSet(merged_df$seq)
names(new_fasta_string) <- merged_df$new_headers
writeXStringSet(new_fasta_string, opt$fastaOutput)



