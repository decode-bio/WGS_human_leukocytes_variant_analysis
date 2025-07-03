#!/bin/bash

# WGS data analysis using the human leucocytes genome
# here I am using the human leucocytes sample from the sra DAtabase

# Step 1: fetching the data from the database with the help of fasterq to reduce one step

fasterq-dump SRR34149094

# Step 2: After donloading the file we will first look for the quality control

fastqc *.fastq

# Step 3: The sequence need not to be trimeed as the quality seems fine. but still for the comparision purpose we\
# still run trimmomatic for the comparison pourpose.
read1=home/sne_desh/rna_seq/P3/SRR34149094_1
read2=home/sne_desh/rna_seq/P3/SRR34149094_2
OP_read1_P=home/sne_desh/rna_seq/P3/SRR34149094_1_P
OP_read1_U=home/sne_desh/rna_seq/P3/SRR34149094_1_U
OP_read2_P=home/sne_desh/rna_seq/P3/SRR34149094_1_P
OP_read2_U=home/sne_desh/rna_seq/P3/SRR34149094_1_U
java -jar $Trimmomatic PE $read1 $read2 \
$OP_read1_P $OP_read1_U $OP_read2_P $OP_read2_U\
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#Step 4: After comparing the both trimmed and untrimmed sequences it is concluded that, the untrimmed sequence is\
# fine to proceed with downstream analysis.
# Now lets run Hisat2 alignement

hisat2 -q -x home/sne_desh/HISAT2/grch38/genome -1 $read1 -2 $read2 -S /sne_desh/WGS/P1/Aligned.sam

# The output showing as terminated, because of the RAM, I have been using 8 GB but it required more. Lets just wait\
# till I get new RAM get installed.

# step 4: Now lets again do fastqc check to see the quality of the trimmed seq

# fastqc *.fastq

# Step 5: Running HISAT2  for the alighment

/usr/bin/hisat2 -q  -x /home/sne_desh/HISAT2/grch38/genome -U /home/sne_desh/rna_seq/P2/SRR34149094.fastq\
-S SRR34149094.sam

