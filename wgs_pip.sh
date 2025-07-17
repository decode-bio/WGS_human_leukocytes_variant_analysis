#!/bin/bash

# WGS data analysis using the human leucocytes sample from the SRA Database

# Step 1: Fetching the data from the database with the help of fasterq (it reduces one step)

fasterq-dump SRR34149094

# Step 2: After downloading the file we will first look for the quality control of sequence.

fastqc *.fastq

# Step 3: The sequence need not to be trimeed as the quality seems fine. But lets compare before and after trimming and then decide.

Trimmomatic=/home/sne_desh/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar
read1=/home/sne_desh/WGS/P1/SRR34149094_1.fastq
read2=/home/sne_desh/WGS/P1/SRR34149094_2.fastq
OP_read1_P=/home/sne_desh/WGS/P1/SRR34149094_1_P.fastq
OP_read1_U=/home/sne_desh/WGS/P1/SRR34149094_1_U.fastq
OP_read2_P=/home/sne_desh/WGS/P1/SRR34149094_1_P.fastq
OP_read2_U=/home/sne_desh/WGS/P1/SRR34149094_1_U.fastq

java -jar $Trimmomatic PE $read1 $read2 $OP_read1_P $OP_read1_U $OP_read2_P $OP_read2_U \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Note: After comparing the both trimmed and untrimmed sequences quality it seems to me that the untrimmed sequence is fine to proceed further.

# Step4: Hisat2 alignement

/usr/bin/hisat2 -q -x /home/sne_desh/HISAT2/grch38/genome -1 /home/sne_desh/WGS/P1/SRR34149094_1.fastq -2 /home/sne_desh/WGS/P1/SRR34149094_2.fastq \
 -S SRR341449094.sam

# The output showing as terminated, because of the RAM, I have been using 8 GB but it required more. Lets just wait\
# till I get new RAM get installed.

# sucessufully run the alignement after expanding the RAM (Finally after two weeks).

# Step 5: Samtool and duplicate removal

# lets convert samfile to bam first.
samtools view -b -S SRR341449094.sam > SRR34149094_unsorted.bam

# lets convert unsorted bam to sorted bam
samtools sort -o SRR34149094_sorted.bam SRR34149094_unsorted.bam

# lets index the sorted bamfiles

samtools index  SRR34149094_sorted.bam

# Statistics of Aligned reads
samtools flagstat SRR34149094_sorted.bam

# Step 6: Variant calling using GATK

