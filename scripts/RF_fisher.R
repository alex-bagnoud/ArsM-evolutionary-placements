#!/usr/bin/env Rscript

# See https://github.com/alex-bagnoud/ArsM-evolutionary-placements/ for details about this script.

library("optparse")

option_list = list(
    make_option(c("-b", "--blastReport"), type="character", default=NULL, 
                help="BLAST report file name, generated using the following output option:  
                '-outfmt 6 sseqid qseqid qstart qend evalue bitscore pident'. The query file of this BLAST
                should be the one specified under to '-td' option.", metavar="character"),
    make_option(c("-t", "--transaltedDB"), type="character", default=NULL, 
                help="Translated database file name. The reading frame should be specified at the end of each
                sequence header, and must be separated from the rest of the header by '_'. No other '_' should
                be present in the sequence header.", metavar="character"),
    make_option(c("-f", "--outFasta"), type="character", default="out.fasta", 
                help="Output fasta file name [default= %default]", metavar="character"),
    make_option(c("-c", "--outBlast"), type="character", default="blast_report.csv", 
                help="Parsed BLAST report output file name (.csv) [default= %default]", metavar="character"),
    make_option(c("-ev", "--Evalue"), type="numeric", default=1e-6, 
                help="BLAST e-value cutoff [default= %default]", metavar="numeric")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Import blast output table file

blast_file <- read.table(opt$blastReport)

## Split the first column into 3

library(stringr)

V2_splitted <- str_split_fixed(blast_file$V2, "_", 2)

## Rename columns

blast_file2 <- data.frame("query_id" = blast_file$V1, "uniq_id" = V2_splitted[,1],
                          "reading_frame" = V2_splitted[,2], "complete_header" = blast_file$V2,
                          "subject_start" = blast_file$V3, "subject_end" = blast_file$V4,
                          "evalue" = blast_file$V5, "score" = blast_file$V6, "p_identity" = blast_file$V7)


## Sort tables and keep the reading frame for each unique sequence

blast_file2 <- blast_file2[order(blast_file2$evalue),]

blast_file3 <- blast_file2[!duplicated(blast_file2$uniq_id),]

## e-value cutoff

blast_file4 <- blast_file3[blast_file3$evalue <= opt$Evalue,]


## Retrieve aa sequences from sequence IDs:

library("Biostrings")

prot_seq = readAAStringSet(opt$transaltedDB)

# split protein header and keep first part only

seq_name_splitted <- strsplit(names(prot_seq), " ")
seq_name <- sapply(seq_name_splitted, `[[`, 1)

sequence = paste(prot_seq)

prot_seq_df <- data.frame(seq_name, sequence)

blast_file5 <- merge(blast_file4, prot_seq_df, by.x = "complete_header", by.y = "seq_name")


## Export files

write.csv(blast_file5, file = opt$outBlast, row.names = FALSE, quote = FALSE)

library(seqinr)

write.fasta(as.list(blast_file5$sequence), as.list(as.character(blast_file5$complete_header)), opt$outFasta, open = "w", nbchar = 200)
