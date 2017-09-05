# ArsM-evolutionary-placements

This script describes the processing of Illumina MiSeq reads from arsM amplicons sequencing. The results were published here:

Matthew C. Reid, Julien Maillard, Alexandre Bagnoud, Leia Falquet, Phu Le Vo, and Rizlan Bernier-Latmani (2017). [**Arsenic Methylation Dynamics in a Rice Paddy Soil Anaerobic Enrichment Culture**](https://dx.doi.org/10.1021/acs.est.7b02970). Environmental Science & Technology. DOI: 10.1021/acs.est.7b02970

#### Table of content

1. [Example](#example)
2. [Example2](#example2)
3. [Third Example](#third-example)

## Example
## Example2
## Third Example



#### Input files

1) Raw fastq files to download from NCBI SRA (see below)

2) ArsM database files, available here: [ArsM-evolutionary-placements/0-input_files/](/0-input_files)


#### Primers used:
330F (5’- GYI WWN GGI VTN GAY ATG A-3’)

385R (5’- ARR TTI AYI ACR CAR TTN S-3’) 


#### Programs used:
* cutadapt v.1.12
* usearch v9.2.64_i86osx32
* blast version 2.2.30+
* EMBOSS v.6.3.1
* R version 3.3.1 (2016-06-21)
* MAFFT version 7: http://mafft.cbrc.jp/alignment/server/
* IQ-TREE 1.5.3: http://iqtree.cibiv.univie.ac.at/
* BMGE v. 1.12 on Galaxy@Pasteur
* BLAST 2.2.30+
* RAxML v.8.2.9
* [RF_fisher.R](scripts/RF_fisher.R) (in-house R script)
* [add_size_to_otus.R](scripts/add_size_to_otus.R) (in-house R script)
* [extract_megan_annotations.R](scripts/extract_megan_annotations.R) (in-house R script)
* seqtk v.1.2 (https://github.com/lh3/seqtk)
* iTOL v.3 (http://itol.embl.de/)
* Inkscape v.0.91
* DIAMOND v.0.8.37
* MEGAN v.6.7.11
* fastq-dump v.2.8.2 from SRA Toolkit

#### 0 Download fastq files from SRA

##### 0.1) Download the fastq file of the paddy soil sample
```
fastq-dump --split-files SRR5972856
mv SRR5972856_1.fastq 0-input_files/SRR5972856_R1.fastq
mv SRR5972856_2.fastq 0-input_files/SRR5972856_R2.fastq
```

##### 0.2) Download the fastq file of the anaerobic enrichment cultures
```
fastq-dump --split-files SRR5972857
mv SRR5972857_1.fastq 0-input_files/SRR5972857_R1.fastq
mv SRR5972857_2.fastq 0-input_files/SRR5972857_R2.fastq
```

#### 1) Define protein OTUs
```
mkdir 1-reid_otus
```

##### 1.1) Merge the paired-end reads
```
usearch -fastq_mergepairs 0-input_files/SRR5972856_R1.fastq -relabel @ -fastqout 1-reid_otus/1-untrimmed_merged.fq
```

##### 1.2) Remove the primer sequences and filter out reads with too many primer mismatches
```
cutadapt -g ^GYNWWNGGNVTNGAYATGA -o 1-reid_otus/2-trimmed_forward.fastq --untrimmed-output 1-reid_otus/2-untrimmed_forward.fastq 1-reid_otus/1-untrimmed_merged.fq

cutadapt -a SNAAYTGYGTNRTNAAYYT$ -o 1-reid_otus/2-trimmed_forward_reverse.fastq --untrimmed-output 1-reid_otus/2-untrimmed_forward_reverse.fastq 1-reid_otus/2-trimmed_forward.fastq
```

##### 1.3) Filter out merged reads with an expect error greater than 1
```
usearch -fastq_filter 1-reid_otus/2-trimmed_forward_reverse.fastq -fastq_maxee 1 -fastaout 1-reid_otus/3-filtered.fasta
```

How many sequences were kept with these quality parameters?
```
for s in 0-input_files/*_R1*.fastq; do
	sample_name=$(echo $s | cut -d "." -f 1 | sed "s/...$//g" | cut -d "/" -f 2);
	count1=$(grep -c "^@${sample_name}" 1-reid_otus/1-untrimmed_merged.fq);	
	count2=$(grep -c "^>${sample_name}" 1-reid_otus/3-filtered.fasta);
	echo -e "$sample_name\t$count1\t$count2";	
done
```
output: ```SRR5972856	171387	166201```


##### 1.4) Dereplication and discargding of singletons
```
usearch -fastx_uniques 1-reid_otus/3-filtered.fasta -fastaout 1-reid_otus/4-1_uniques.fasta -relabel Uniq -sizeout
usearch -sortbysize 1-reid_otus/4-1_uniques.fasta -fastaout 1-reid_otus/4-2_uniques_nosingle.fasta -minsize 2
```

##### 1.5) Translation of uniques sequences to proteins sequences
```
transeq 1-reid_otus/4-2_uniques_nosingle.fasta 1-reid_otus/5_uniques_translated_6frames.fasta -table 11 -frame=6
```

##### 1.6) Fishing the correct reading frames by blasting translated uniques reads on arSM reference sequence
```
cat 0-input_files/arsm_11protdb_char.fasta 0-input_files/arsm_723protdb_non-char.fasta > 1-reid_otus/6-1_arsm_734protdb.fasta
makeblastdb -in 1-reid_otus/6-1_arsm_734protdb.fasta -out 1-reid_otus/6-2_arsM_blast_db -dbtype prot
blastp -query 1-reid_otus/5_uniques_translated_6frames.fasta -db 1-reid_otus/6-2_ArsM_blast_db -out 1-reid_otus/6-3_blast_report.tab -outfmt "6 sseqid qseqid qstart qend evalue bitscore pident" -max_target_seqs 1
```

##### 1.7) BLAST output parsing
This in-house R script parses to previous BLAST ouput, selects for each translated sequence which of the reading frame matches ArsM proteins, and extract the correct proteins sequences into a fasta file.
```
scripts/RF_fisher.R -b 1-reid_otus/6-3_blast_report.tab -t 1-reid_otus/5_uniques_translated_6frames.fasta -f 1-reid_otus/7-prot.fasta -c 1-reid_otus/7-parsed_blast_report.csv
```

##### 1.8) Extract uniques proteins sequences
```
usearch -fastx_uniques 1-reid_otus/7-prot.fasta -fastaout 1-reid_otus/8-prot_uniques.fasta -relabel uniq_prot -sizeout
```

##### 1.9) Define clusters with 97%-similarity threshold
```
usearch -cluster_otus 1-reid_otus/8-prot_uniques.fasta -otus 1-reid_otus/9-prot_otus.fasta -relabel Otu
```

##### 1.10) Translation of the reads to protein sequences (in all 6 reading frames)
```
transeq 1-reid_otus/2-trimmed_forward_reverse.fastq 1-reid_otus/10-trimmed_forward_reverse_translated.fasta -table 11 -frame=6
```

##### 1.11) Build OTU table (based on 97%-similarity threshold)
```
usearch -usearch_global 1-reid_otus/10-trimmed_forward_reverse_translated.fasta -db 1-reid_otus/9-prot_otus.fasta -id .97 -otutabout 1-reid_otus/11-prot_otu_table.txt
```

##### 1.12) Add size information for each OTU
```
scripts/add_size_to_otus.R -t 1-reid_otus/11-prot_otu_table.txt -i 1-reid_otus/9-prot_otus.fasta -f 1-reid_otus/12-prot_otus_size.fasta
```

##### 1.13) Cluster protein OTUs sequences at 90% similarity
```
usearch -cluster_fast 1-reid_otus/12-prot_otus_size.fasta -id 0.9 -sort size -centroids 1-reid_otus/13-prot_otus_uc90.fasta
```


#### 2) Compiling arsM protein OTUs database from Jia's primers
```
mkdir 2-jia_otus
```

##### 2.0) Recover arsM OTUs from Jia's primers:

Study 1:
  * Jia et al., 2013, ES&T, Microbial Arsenic Methylation in Soil and Rice Rhizosphere
  * NCBI accession numbers: JQ924198 to JQ924282
  * Sequenced retrieved using NCBI Batch Entrez

Study 2:
  * Zhang et al., 2015, ES&T, Diversity and Abundance of Arsenic Biotransformation Genes in Paddy Soils from Southern China
  * NCBI accession numbers: KP060920 to KP060962
  * Sequenced retrieved using NCBI Batch Entrez

Study 3:
  * Desoeuvre et al., 2016, Microbial Ecology, Diversity and Distribution of Arsenic-Related Genes Along a Pollution Gradient in a River Affected by Acid Mine Drainage
  * NCBI accession numbers: KR051692 to KR051724
  * Sequenced retrieved using NCBI Batch Entrez

Sequences saved as: ```0-input_files/jia_otus.fasta```


##### 2.1) translate dna otus into prot otus
```
transeq 0-input_files/jia_otus.fasta  2-jia_otus/1-jia_otus_translated_6frames.fasta -table 11 -frame=6
```

##### 2.2) Filter out correct reading frame
```
blastp -query 2-jia_otus/1-jia_otus_translated_6frames.fasta -db 1-reid_otus/6-2_ArsM_blast_db -out 2-jia_otus/2-zia_prot_otus_fished_by_arsMdb.txt -outfmt "6 sseqid qseqid qstart qend evalue bitscore pident" -max_target_seqs 1
```

##### 2.3) BLAST output parsing
```
scripts/RF_fisher.R -b 2-jia_otus/2-zia_prot_otus_fished_by_arsMdb.txt -t 2-jia_otus/1-jia_otus_translated_6frames.fasta -f 2-jia_otus/3-jia_prot.fasta -c 2-jia_otus/3-jia_parsed_blast_report.csv
```

##### 2.4) Cluster Jia OTUs sequences at 90% similarity
```
usearch -cluster_fast 2-jia_otus/3-jia_prot.fasta -id 0.9 -sort length -centroids 2-jia_otus/4-jia_prot_uc90.fasta
```



#### 3) Computing a arsM phylogenetic tree
```
mkdir 3-arsm_ref_tree
```

##### 3.1) Cat both protein OTUs files
```
cat 1-reid_otus/13-prot_otus_uc90.fasta 2-jia_otus/4-jia_prot_uc90.fasta > 3-arsm_ref_tree/1-arsm_prot_otus_uc90.fasta
```

##### 3.2) Re-order databases sequences
```
usearch -sortbylength 0-input_files/arsm_11protdb_char.fasta -fastaout 3-arsm_ref_tree/2-1_arsm_11protdb_sorted.fasta
sed 's/>.*/&;size=2;/' 3-arsm_ref_tree/2-1_arsm_11protdb_sorted.fasta > 3-arsm_ref_tree/2-2_arsm_11protdb_sorted_labeled.fasta
usearch -sortbylength 0-input_files/arsm_723protdb_non-char.fasta -fastaout 3-arsm_ref_tree/2-3_arsm_723protdb_sorted.fasta
```

##### 3.3) Center the ArsM db of uncharacterized proteins using a similarity treshold of 80% with ArsM OTUs:
```
usearch -usearch_global 3-arsm_ref_tree/2-3_arsm_723protdb_sorted.fasta -db 3-arsm_ref_tree/1-arsm_prot_otus_uc90.fasta -id 0.8 -matched 3-arsm_ref_tree/3-centered80_arsm_723protdb.fasta
```

##### 3.4) Add biochemically-characterized sequences to the new set of reference seqences
```
cat 3-arsm_ref_tree/2-2_arsm_11protdb_sorted_labeled.fasta 3-arsm_ref_tree/3-centered80_arsm_723protdb.fasta > 3-arsm_ref_tree/4-centered80_arsm_734protdb.fasta
```


##### 3.5) Cluster centered database at 90% similarity
```
usearch -cluster_fast 3-arsm_ref_tree/4-centered80_arsm_734protdb.fasta -id 0.9 -sort size -centroids 3-arsm_ref_tree/5-centered80_arsm_protdb_uc90.fasta
```

##### 3.6) Alignement with mafft
* http://mafft.cbrc.jp/alignment/software/
* MAFFT version 7
* Default parameters
* input: ```3-arsm_ref_tree/5-centered80_arsm_protdb_uc90.fasta```
* output: ```3-arsm_ref_tree/6-centered80_arsm_protdb_uc90_mafft.fasta```


##### 3.7) Alignement parsing with BMGE
* Galaxy@Pasteur (https://galaxy.pasteur.fr/)
* BMGE v. 1.12
* Sequence Coding	AA	
* matrice	BLOSUM	
* Estimated matrice BLOSUM	62	
* Sliding Windows Size	3	
* Maximum entropy threshold	1.0	
* Gap Rate cut-off [ 0-1 ]	0.9	
* Minimum Block Size	1
* Output file(s)	Fasta
* input: ```3-arsm_ref_tree/6-centered80_arsm_protdb_uc90_mafft.fasta```
* output: ```3-arsm_ref_tree/7-centered80_arsm_protdb_uc90_bmge.fasta```


##### 3.8) IQ-TREE
* http://iqtree.cibiv.univie.ac.at/
* version 1.5.3
* default parameters
* input: ```3-arsm_ref_tree/7-centered80_arsm_protdb_uc90_bmge.fasta```
* outputs: ```3-arsm_ref_tree/8*```




#### 4) EPA analysis of Jia's OTUs
```
mkdir 4-epa_jia_otus
```

##### 4.1) Alignement with mafft
* http://mafft.cbrc.jp/alignment/software/
* MAFFT version 7
* ```--addfragments``` option with default parameters
* input file 1 (bmge-parsed alignement): ```3-arsm_ref_tree/7-centered80_arsm_protdb_uc90_bmge.fasta```
* input file 2 (otus sequences): ```2-jia_otus/4-jia_prot_uc90.fasta```
* output: ```4-epa_jia_otus/1-1-jia_otus_added_to_ref_alignment_mafft.fasta```

Parse the file to make it RAxML-compatible
```
less 4-epa_jia_otus/1-1-jia_otus_added_to_ref_alignment_mafft.fasta | tr '[|]:)("/,' '_' > 4-epa_jia_otus/1-2-jia_otus_added_to_ref_alignment_mafft.fasta
```

##### 4.2) EPA analysis with RAxML
Be sure the select the  substitution model that was selected by IQ-TREE. In this case, LG+F+I+G4 IQ-TREE model corresponds to PROTGAMMAILGF RAxML model.
```
raxmlHPC -f v -s 4-epa_jia_otus/1-2-jia_otus_added_to_ref_alignment_mafft.fasta -t 3-arsm_ref_tree/8-centered80_arsm_protdb_uc90_bmge.fasta.contree -m PROTGAMMAILGF --epa-keep-placements 1 -n epa_jia
```

Move manually and rename all output files	this way:  ```4-epa_jia_otus/2-*```



#### 5) EPA analysis of Reid's OTUs
```
mkdir 5-epa_reid_otus
```

##### 5.1) Alignement with mafft
* http://mafft.cbrc.jp/alignment/software/
* MAFFT version 7
* ```--addfragments``` option with default parameters
* input file 1 (bmge-parsed alignement): ```3-arsm_ref_tree/7-centered80_arsm_protdb_uc90_bmge.fasta```
* input file 2 (otus sequences): ```1-reid_otus/13-prot_otus_uc90.fasta```
* output: ```5-epa_reid_otus/1-1-reid_otus_added_to_ref_alignment_mafft.fasta```

Parse the file to make it RAxML-compatible
```
less 5-epa_reid_otus/1-1-reid_otus_added_to_ref_alignment_mafft.fasta | tr '[|]:)("/,;' '_' > 5-epa_reid_otus/1-2-reid_otus_added_to_ref_alignment_mafft.fasta
```

##### 5.2) EPA analysis with RAxML
Be sure the select the  substitution model that was selected by IQ-TREE. In this case, LG+F+I+G4 IQ-TREE model corresponds to PROTGAMMAILGF RAxML model.
```
raxmlHPC -f v -s 5-epa_reid_otus/1-2-reid_otus_added_to_ref_alignment_mafft.fasta -t 3-arsm_ref_tree/8-centered80_arsm_protdb_uc90_bmge.fasta.contree -m PROTGAMMAILGF --epa-keep-placements 1 -n epa_reid
```

Move manually and rename all output files	this way:  ```5-reid_jia_otus/2-*```



#### 6) Downscaling Reid's reads
Because of the huge discrepancy between the number of reads from Reid's sequencing (166'197 trimmed reads) and the reads from other sequencing projects based on Jia's primers sets (531 clones), 531 Reid's reads were randomly selected, before proceeding with the rest of the pipeline. This intend to increase the comparability of the diversity spectrum of the 2 primers sets.

```
mkdir 6-downscaling_test
```


##### 6.1) Use seqtk (v.1.2) to subsample the trimmed raw reads
```
seqtk sample -s100 1-reid_otus/3-filtered.fasta 531 > 6-downscaling_test/1-531subset_reads.fasta
```

##### 6.2) Dereplication and discargding of singletons
```
usearch -fastx_uniques 6-downscaling_test/1-531subset_reads.fasta -fastaout 6-downscaling_test/2-uniques.fasta -relabel Uniq -sizeout
```

##### 6.3) Translation of uniques sequences to proteins sequences
```
transeq 6-downscaling_test/2-uniques.fasta 6-downscaling_test/3-uniques_translated_6frames.fasta -table 11 -frame=6
```

##### 6.4) Fishing the correct reading frames by blasting translated uniques reads on arSM reference sequence
```
blastp -query 6-downscaling_test/3-uniques_translated_6frames.fasta -db 1-reid_otus/6-2_ArsM_blast_db -out 6-downscaling_test/4-blast_report.tab -outfmt "6 sseqid qseqid qstart qend evalue bitscore pident" -max_target_seqs 1
```

##### 6.5) BLAST output parsing
This in-house R script parses to previous BLAST ouput, selects for each translated sequence which of the reading frame matches ArsM proteins, and extract the correct proteins sequences into a fasta file.
```
scripts/RF_fisher.R -b 6-downscaling_test/4-blast_report.tab -t 6-downscaling_test/3-uniques_translated_6frames.fasta -f 6-downscaling_test/5-prot.fasta -c 6-downscaling_test/5-parsed_blast_report.csv
```
##### 6.6) Extract uniques proteins sequences
```
usearch -fastx_uniques 6-downscaling_test/5-prot.fasta -fastaout 6-downscaling_test/6-prot_uniques.fasta -relabel uniq_prot -sizeout
```

##### 6.7) Define clusters with 97%-similarity threshold
```
usearch -cluster_otus 6-downscaling_test/6-prot_uniques.fasta -otus 6-downscaling_test/7-prot_otus.fasta -relabel Otu
```

##### 6.8) Translation of the reads to protein sequences (in all 6 reading frames)
```
transeq 6-downscaling_test/1-531subset_reads.fasta 6-downscaling_test/8-t531subset_translated.fasta -table 11 -frame=6
```

##### 6.9) Build OTU table (based on 97%-similarity threshold)
```
usearch -usearch_global 6-downscaling_test/8-t531subset_translated.fasta -db 6-downscaling_test/7-prot_otus.fasta -id .97 -otutabout 6-downscaling_test/9-prot_otu_table.txt
```

##### 6.10) Add size information for each OTU
```
scripts/add_size_to_otus.R -t 6-downscaling_test/9-prot_otu_table.txt -i 6-downscaling_test/7-prot_otus.fasta -f 6-downscaling_test/10-prot_otus_size.fasta
```

##### 6.11) Cluster protein OTUs sequences at 90% similarity
```
usearch -cluster_fast 6-downscaling_test/10-prot_otus_size.fasta -id 0.9 -sort size -centroids 6-downscaling_test/11-prot_otus_uc90.fasta
```

##### 6.12) Alignement with mafft
* http://mafft.cbrc.jp/alignment/software/
* MAFFT version 7
* ```--addfragments``` option with default parameters
* input file 1 (bmge-parsed alignement): ```3-arsm_ref_tree/7-centered80_arsm_protdb_uc90_bmge.fasta```
* input file 2 (otus sequences): ```6-downscaling_test/11-prot_otus_uc90.fasta```
* output: ```6-downscaling_test/12-addfragments_mafft.fasta```

##### 6.13) Parse the file to make it RAxML-compatible
```
less 6-downscaling_test/12-addfragments_mafft.fasta | tr '[|]:)("/,;' '_' > 6-downscaling_test/13-addfragments_mafft.fasta
```

##### 6.14) EPA analysis with RAxML
Be sure the select the  substitution model that was selected by IQ-TREE. In this case, LG+F+I+G4 IQ-TREE model corresponds to PROTGAMMAILGF RAxML model.
```
raxmlHPC -f v -s 6-downscaling_test/13-addfragments_mafft.fasta -t 3-arsm_ref_tree/8-centered80_arsm_protdb_uc90_bmge.fasta.contree -m PROTGAMMAILGF --epa-keep-placements 1 -n epa_reid_ds
```

Move manually and rename all output files	this way:  ```6-downscaling_test/14-*```

#### 7) Tree annotation

Trees were manually annotated with the interactive tree of life (iTol v.3) on http://itol.embl.de/. They were then merged and further annotated with Inkscape v.091

```
mkdir 7-itol_annotations
```



#### 8) Define protein OTUs
```
mkdir 8-7samples_otus
```

##### 8.1) Merge the paired-end reads
```
usearch -fastq_mergepairs 0-input_files/*_R1*.fastq -relabel @ -fastqout 8-7samples_otus/1-untrimmed_merged.fq
```

##### 8.2) Remove the primer sequences and filter out reads with too many primer mismatches
```
cutadapt -g ^GYNWWNGGNVTNGAYATGA -o 8-7samples_otus/2-trimmed_forward.fastq --untrimmed-output 8-7samples_otus/2-untrimmed_forward.fastq 8-7samples_otus/1-untrimmed_merged.fq

cutadapt -a SNAAYTGYGTNRTNAAYYT$ -o 8-7samples_otus/2-trimmed_forward_reverse.fastq --untrimmed-output 8-7samples_otus/2-untrimmed_forward_reverse.fastq 8-7samples_otus/2-trimmed_forward.fastq
```

##### 8.3) Filter out merged reads with an expect error greater than 1
```
usearch -fastq_filter 8-7samples_otus/2-trimmed_forward_reverse.fastq -fastq_maxee 1 -fastaout 8-7samples_otus/3-filtered.fasta
```

How many sequences were kept with these quality parameters?
```
for s in 0-input_files/*_R1*.fastq; do
	sample_name=$(echo $s | cut -d "." -f 1 | sed "s/...$//g" | cut -d "/" -f 2);
	count1=$(grep -c "^@${sample_name}" 8-7samples_otus/1-untrimmed_merged.fq);	
	count2=$(grep -c "^>${sample_name}" 8-7samples_otus/3-filtered.fasta);
	echo -e "$sample_name\t$count1\t$count2";	
done
```
output: ```SRR5972856	171387	166201```
```SRR5972857	342190	328212```


##### 8.4) Dereplication and discargding of singletons
```
usearch -fastx_uniques 8-7samples_otus/3-filtered.fasta -fastaout 8-7samples_otus/4-1_uniques.fasta -relabel Uniq -sizeout
usearch -sortbysize 8-7samples_otus/4-1_uniques.fasta -fastaout 8-7samples_otus/4-2_uniques_nosingle.fasta -minsize 2
```

##### 8.5) Translation of uniques sequences to proteins sequences
```
transeq 8-7samples_otus/4-2_uniques_nosingle.fasta 8-7samples_otus/5_uniques_translated_6frames.fasta -table 11 -frame=6
```

##### 8.6) Fishing the correct reading frames by blasting translated uniques reads on arSM reference sequence
```
cat 0-input_files/arsm_11protdb_char.fasta 0-input_files/arsm_723protdb_non-char.fasta > 8-7samples_otus/6-1_arsm_734protdb.fasta
makeblastdb -in 8-7samples_otus/6-1_arsm_734protdb.fasta -out 8-7samples_otus/6-2_arsM_blast_db -dbtype prot
blastp -query 8-7samples_otus/5_uniques_translated_6frames.fasta -db 8-7samples_otus/6-2_ArsM_blast_db -out 8-7samples_otus/6-3_blast_report.tab -outfmt "6 sseqid qseqid qstart qend evalue bitscore pident" -max_target_seqs 1
```

##### 8.7) BLAST output parsing
This in-house R script parses to previous BLAST ouput, selects for each translated sequence which of the reading frame matches ArsM proteins, and extract the correct proteins sequences into a fasta file.
```
scripts/RF_fisher.R -b 8-7samples_otus/6-3_blast_report.tab -t 8-7samples_otus/5_uniques_translated_6frames.fasta -f 8-7samples_otus/7-prot.fasta -c 8-7samples_otus/7-parsed_blast_report.csv
```

##### 8.8) Extract uniques proteins sequences
```
usearch -fastx_uniques 8-7samples_otus/7-prot.fasta -fastaout 8-7samples_otus/8-prot_uniques.fasta -relabel uniq_prot -sizeout
```

##### 8.9) Define clusters with 97%-similarity threshold
```
usearch -cluster_otus 8-7samples_otus/8-prot_uniques.fasta -otus 8-7samples_otus/9-prot_otus.fasta -relabel Otu
```

##### 8.10) Translation of the reads to protein sequences (in all 6 reading frames)
```
transeq 8-7samples_otus/2-trimmed_forward_reverse.fastq 8-7samples_otus/10-trimmed_forward_reverse_translated.fasta -table 11 -frame=6
```

##### 8.11) Build OTU table (based on 97%-similarity threshold)
```
usearch -usearch_global 8-7samples_otus/10-trimmed_forward_reverse_translated.fasta -db 8-7samples_otus/9-prot_otus.fasta -id .97 -otutabout 8-7samples_otus/11-prot_otu_table.txt
```


##### 8.12) Create DIAMOND database

Merge 2 database fasta files
```
cat 0-input_files/arsm_723protdb_non-char.fasta 0-input_files/arsm_11protdb_char.fasta | tr '_bacterium' '_' > 0-input_files/arsm_734prot_db.fasta
```
Remove manually in ```0-input_files/arsm_734prot_db.fasta``` the two endosymbionts and save it as ```0-input_files/arsm_732prot_db.fasta````

Create a diamond database
```
diamond makedb --in 0-input_files/arsm_732prot_db.fasta --db 8-7samples_otus/12-arsm_732prot_db
```


##### 8.13) Blast proteins OTUs on arsm protein database
```
diamond blastp --query 8-7samples_otus/9-prot_otus.fasta --db 8-7samples_otus/12-arsm_732prot_db --daa 8-7samples_otus/13-1-reid_otus_diamond
diamond view -a 8-7samples_otus/13-1-reid_otus_diamond.daa > 8-7samples_otus/13-2-reid_otus_diamond.txt
```


##### 8.14) LCA algorithm with MEGAN

* Open MEGAN 6
* File > Meganize DAA File
* Files: ```8-tax_annotation/13-1-reid_otus_diamond```
* Taxonomy: Tick 'Parse Taxon Names' and 'Use Id parsing'
* LCA algorithm: Weighted algorithm with default paramaters
* Open meganized file
* Select all
* Click on 'Inspect the read-to-functions assignment'
* Expand all taxons
* Copy paste, and save as ```8-tax_annotation/14-otu_annotations_megan.txt```

##### 8.7) Analysis the MEGAN output file

Run ```extract_megan_annotations.R``` on R, whose output is an OTU table that contains taxonomic annotations, saved as ```8-7samples_otus/15-prot_otu_table_tax.txt```.


The full taxonomic path of OTUs annotated as 'Actinobacteria <phylum>' by LCA could not be retrieved. For those, the taxonomic path was manually written in ```8-7samples_otus/15-prot_otu_table_tax.txt```. It concerned Otu2750, Otu3624, Otu4069, Otu4285, Otu4713, Otu5012, Otu5535, Otu604, Otu6046, Otu6563, Otu7358, Otu7558, Otu8097, and Otu8763














