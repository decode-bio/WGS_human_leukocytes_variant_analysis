#!/bin/bash

# WGS data analysis using the human leucocytes sample from the SRA Database

# Step 1: Fetching the data from the database with the help of fasterq (it reduces one step)

#fasterq-dump SRR34149094

# Step 2: After downloading the file we will first look for the quality control of sequence.

#fastqc *.fastq

# Step 3: The sequence need not to be trimeed as the quality seems fine. But lets compare before and after trimming and then decide.

#Trimmomatic=/home/sne_desh/Tools/Trimmomatic-0.39/trimmomatic-0.39.jar
#read1=/home/sne_desh/WGS/P1/SRR34149094_1.fastq
#read2=/home/sne_desh/WGS/P1/SRR34149094_2.fastq
#OP_read1_P=/home/sne_desh/WGS/P1/SRR34149094_1_P.fastq
#OP_read1_U=/home/sne_desh/WGS/P1/SRR34149094_1_U.fastq
#OP_read2_P=/home/sne_desh/WGS/P1/SRR34149094_1_P.fastq
#OP_read2_U=/home/sne_desh/WGS/P1/SRR34149094_1_U.fastq

#java -jar $Trimmomatic PE $read1 $read2 $OP_read1_P $OP_read1_U $OP_read2_P #$OP_read2_U \
#ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Note: After comparing the both trimmed and untrimmed sequences quality it seems to me that the untrimmed sequence is fine to proceed further.

# Step4: As I am doing the WGS, so we use the BWA alignement

#After installing bwa and downloading the reference genome and indexing it.

#bwa mem \
#/home/sne_desh/Ref/Homo_sapiens_assembly38.fasta \
#/home/sne_desh/WGS/P1/SRR34149094_1.fastq \
#/home/sne_desh/WGS/P1/SRR34149094_2.fastq \
#> aligned.sam


# Step 5: Samtool and duplicate removal

# lets convert samfile to bam first.
# samtools view -b -S aligned.sam > aligned_unsorted.bam

# lets convert unsorted bam to sorted bam
# samtools sort -n -o aligned_sorted.bam aligned_unsorted.bam

#Fix mate-pair information

# samtools fixmate -m aligned_sorted.bam aligned_sorted.fixmate.bam


# sort by genomic coordinates
# samtools sort -o aligned_position.fixmate.bam aligned_sorted.fixmate.bam

# mark and remove duplicates
# samtools markdup -r aligned_position.fixmate.bam aligned_sorted.rm.bam

# lets index the final sorted bamfiles

# samtools index  aligned_sorted.rm.bam

# Statistics of Aligned reads
# samtools flagstat aligned_sorted.rm.bam

# Step 6: Variant calling using GATK
#ref=/home/sne_desh/Ref/Homo_sapiens_grch38.fa
#python3 ~/Tools/gatk-4.6.2.0/gatk \
  # --java-options "-Xmx4g" HaplotypeCaller  \
  # -R $ref \
   #-I aligned_sorted.rm.bam\
   #-O SRR34149094.g.vcf.gz \
  # -ERC GVCF

