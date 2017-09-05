#!/usr/bin/env Rscript

# This script help to extract for each read the full taxonomic path of MEGAN LCA annotations
# The input file must be created by copying-pasting the Inspector windows, after uncollapsing each taxonomic node.
# This script uses a function than retrieve the full taxonomic path from NCBI from any taxonomic node.
# Sometimes, a conflict occur, and one has to select manually the taxonomic path.
# See https://github.com/alex-bagnoud/ArsM-evolutionary-placements/ for more details about this script.


# Set variables

megan_file <- "8-7samples_otus/14-otu_annotations_megan.txt"
label <- "Otu"
otu_table <- "8-7samples_otus/11-prot_otu_table.txt"
blast_file <- "8-7samples_otus/13-2-reid_otus_diamond.txt"
output <- "8-7samples_otus/15-prot_otu_table_tax.txt"

# Import the MEGAN file as two list of the same length,
# one that has the seqeunce header, and another one that has the annotations
con  <- file(megan_file, open = "r")

tax_list <- character()
read_list <- character()

while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if ( !grepl(label, oneLine) ) {
        tax <- oneLine
    }
    else {
        tax_list <- c(tax_list, tax)
        read_list <- c(read_list, oneLine)
    }
} 

close(con)

# Remove brackets from vectors elements
library("stringr")
tax_list <- str_split_fixed(tax_list, " \\[", 2)[,1]
read_list <- str_split_fixed(read_list, " \\[", 2)[,1]

# Merge these 2 vectors into a dataframe

read_tax_df <- data.frame("otu" = read_list, "tax" = tax_list)

# Retriev full taxonomic path from NCBI
library("myTAI")

fullTaxRank <- function(org_name){
    
    tax <- taxonomy(org_name,db = "ncbi")
    
    output <- list("superkingdom" = NA, "phylum" = NA, "class" = NA, "order" = NA, "family" = NA, "genus" = NA, "species" = NA)
    
    tax.d <- tax[tax$rank == "superkingdom",1]
    tax.p <- tax[tax$rank == "phylum",1]
    tax.c <- tax[tax$rank == "class",1]
    tax.o <- tax[tax$rank == "order",1]
    tax.f <- tax[tax$rank == "family",1]
    tax.g <- tax[tax$rank == "genus",1]
    tax.s <- tax[tax$rank == "species",1]
    
    if (!(length(tax.s) == 0)) {
        output["species"] <- tax.s
    }
    if (!(length(tax.g) == 0)) {
        output["genus"] <- tax.g
    }
    if (!(length(tax.f) == 0)) {
        output["family"] <- tax.f
    }
    if (!(length(tax.o) == 0)) {
        output["order"] <- tax.o
    }
    if (!(length(tax.c) == 0)) {
        output["class"] <- tax.c
    }
    if (!(length(tax.p) == 0)) {
        output["phylum"] <- tax.p
    }
    if (!(length(tax.d) == 0)) {
        output["superkingdom"] <- tax.d
    }
    return(unlist(output))
}

full_tax_list <- character()
full_tax_path <- character()
l.1 <- character()
l.2 <- character()
l.3 <- character()
l.4 <- character()
l.5 <- character()
l.6 <- character()
l.7 <- character()

for (i in read_tax_df[,2]) {
    if ( !(i %in% full_tax_list) ) {
        full_tax_list <- c(full_tax_list, i)
        full_tax_path <- fullTaxRank(i)
        l.1 <- c(l.1, full_tax_path[1])
        l.2 <- c(l.2, full_tax_path[2])
        l.3 <- c(l.3, full_tax_path[3])
        l.4 <- c(l.4, full_tax_path[4])
        l.5 <- c(l.5, full_tax_path[5])
        l.6 <- c(l.6, full_tax_path[6])
        l.7 <- c(l.7, full_tax_path[7])
    }
    else {
        l.1 <- c(l.1, full_tax_path[1])
        l.2 <- c(l.2, full_tax_path[2])
        l.3 <- c(l.3, full_tax_path[3])
        l.4 <- c(l.4, full_tax_path[4])
        l.5 <- c(l.5, full_tax_path[5])
        l.6 <- c(l.6, full_tax_path[6])
        l.7 <- c(l.7, full_tax_path[7])
    }
}

read_tax_df$lca_superkingdom <- l.1
read_tax_df$lca_phylum <- l.2
read_tax_df$lca_class <- l.3
read_tax_df$lca_order <- l.4
read_tax_df$lca_family <- l.5
read_tax_df$lca_genus <- l.6


# Import DIAMOND blast table
blast <- read.table(blast_file, header = FALSE)

# Sort dataframe by increasing e-value and decreasing pident
sorted_blast <- blast[order(blast$V11, -blast$V3),]

# Keep the best hit for each OTU
best_hit <- sorted_blast[!duplicated(sorted_blast$V1),]

# Merge annotation and blast dataframes
merged_df <- merge.data.frame(read_tax_df, best_hit, by.x = "otu", by.y = "V1", all = TRUE)

merged_df2 <- merged_df[,c(1,3,4,5,6,7,8,9,10,18)]
names(merged_df2)[c(8,9,10)] <- c("diamond_best_hit", "diamond_pident", "diamond_e-value")

# Import OTU table a merged it with the taxonomix assignemetn dataframe

otu_tab <- read.table(otu_table, header = FALSE, sep = '\t')
names(otu_tab) <- c("OTU_id", "S1", "S2", "S3", "S4", "S5", "S6", "S7")

otu_tab_tax <- merge(otu_tab, merged_df2, by.x = "OTU_id", by.y = "otu", all = TRUE)

# Save merged dataframe as file

write.table(otu_tab_tax, file = output, quote = FALSE, sep = '\t', row.names = FALSE)

