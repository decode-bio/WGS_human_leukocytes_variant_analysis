# Objetive

# Tools
* SRA Toolkit
* Trimmomatic
* Hisat2
* SAMTOOLS
* GATK
* ANNOVAR

# About Data

* Source: GENOMIC
* Organism: Homo sapiens
* Name: MycoTitrate2wgsvXR28q1_S65 WGS seq 1
* Instrument: Illumina NovaSeq 6000
* Strategy: WGS
* Layout: PAIRED

![image](https://github.com/user-attachments/assets/df0ff94a-910a-4596-b33d-f5e4f8f443e8)

# Step1: Downloading the data from SRA database
```
fasterq -dump SRR34149094
```
![image](https://github.com/user-attachments/assets/096275da-033e-478e-8352-f9e186d1b9f2)

# Step2: Quality Control
```
fastqc *.fastq

```
![image](https://github.com/user-attachments/assets/f7b4d773-b008-4f01-92e9-b85fdb45be6a)

# Step3: Quality control analysis
![image](https://github.com/user-attachments/assets/ae767346-da50-4b5b-acc5-f1d3201fe360)
![image](https://github.com/user-attachments/assets/1964bf56-49da-4612-8cde-e81d871d840a)
As the data does not contain the adapter contamination, also the alignment (HISAT2) we re using have smooth trimming. we do not need to perform the trimming of the data as it will reduce our sequence length futher, which we donot want at this point. Hence lets proceed with untrimmed data.
Still for comparision lets run Trimmomatic and compare the trimmed and untrimmed results to get a clear idea.
* Trimmomatic
```
Trimmomatic=/usr/share/java/trimmomatic-0.39.jar
read1=/home/sne_desh/rna_seq/P3/SRR34149094_1.fastq
read2=/home/sne_desh/rna_seq/P3/SRR34149094_2.fastq
OP_read1_P=/home/sne_desh/rna_seq/P3/SRR34149094_1_paired.fastq
OP_read1_U=/home/sne_desh/rna_seq/P3/SRR34149094_1_unpaired.fastq
OP_read2_P=/home/sne_desh/rna_seq/P3/SRR34149094_2_paired.fastq
OP_read2_U=/home/sne_desh/rna_seq/P3/SRR34149094_2_unpaired.fastq

java -jar $Trimmomatic PE $read1 $read2 \
$OP_read1_P $OP_read1_U $OP_read2_P $OP_read2_U\
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
![image](https://github.com/user-attachments/assets/22c0029b-2b1c-4350-ab11-7fbcf024b5b7)




