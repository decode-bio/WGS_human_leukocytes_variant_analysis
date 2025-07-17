# Objetive

# Tools
* SRA Toolkit
* Trimmomatic
* Hisat2
* SAMTOOLS
* SUBREAD
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

java -jar $Trimmomatic PE $read1 $read2 \
$OP_read1_P $OP_read1_U $OP_read2_P $OP_read2_U\
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```
![image](https://github.com/user-attachments/assets/22c0029b-2b1c-4350-ab11-7fbcf024b5b7)

* QC report after trimming
![image](https://github.com/user-attachments/assets/11571411-ca27-4892-92ad-69c1c4f3f57e)
![image](https://github.com/user-attachments/assets/35af96a5-57d2-4cc3-8f31-a0b52a92b250)

* On comparing both the reports, it is clearly seen that the sequence base lenght has been reduced, which will increase the false positive alignment. So I will proceed further with the untrimmed fastq files for the downstream analysis. Lets align our data to the genome using HISAT2.

# Step4: Alignment using HISAT2
Now after QC, we proceed with our original reads and run alignment with HISAT2. For this step the HISAT2 and genome index has already been downloaded.
```
/usr/bin/hisat2 -q -x /home/sne_desh/HISAT2/grch38/genome -1 /home/sne_desh/WGS/P1/SRR34149094_1.fastq -2 /home/sne_desh/WGS/P1/SRR34149094_2.fastq -S SRR34149094.sam
```
<img width="1815" height="358" alt="image" src="https://github.com/user-attachments/assets/e5e5d2b9-50d0-40b0-af5c-009ffb35fb07" />
Our Overall alignement rate is 95.46% which is a good alignment rate and it suggest that our library is free from contamination or adapter issue. Now we get our aligned read in .sam format. Lets continue with the next step to do sorting and duplicate removal.

# Step 5: Sorting and duplicate removal using samtools
* Why Sort BAM Files?
- Sorting aligns reads by genomic coordinates, which is essential for indexing and efficient retrieval.
- Many downstream tools—such as variant callers (GATK), quantifiers (featureCounts), and genome browsers (IGV) require coordinate-sorted BAM files to function correctly.
- A sorted BAM reduces memory overhead during region-based queries, since alignments for any locus are grouped together on disk.
* What Duplicate Removal Does?
- Duplicate removal filters out PCR or optical duplicates—reads that originate from the same original DNA fragment but were sequenced multiple times.
- Leaving duplicates in your data can inflate read counts at certain loci, biasing expression estimates or variant allele frequencies.

Conversion from .sam to .bam
```
samtools -b -S SRR34149094.sam > SRR34149094.bam
```
Sorting the unsorted .bam file
```
samtools sort -o SRR34149094_sorted.bam SRR34149094.bam
```
Indexing the sorted .bam file
```
samtools index SRR34149094_sorted.bam
```
Output
<img width="1778" height="151" alt="image" src="https://github.com/user-attachments/assets/33b518bc-9241-4c97-8d0f-3cdbf6be781c" />
Now I am done with the sorting and removal of duplicate and indexing of my bam file. lets go to the next step i.e. feature count.

# step 6: Feature counts
As for performing the feature counts we will be needing the SUBREAD, so first lets install the feature count to our sysytem by writing the simple command.
```
sudo apt install subread
```
OUTPUT
<img width="1180" height="416" alt="image" src="https://github.com/user-attachments/assets/69ee5911-b0df-4131-9235-031964bdd4d9" />

Now lets perform the feature counts

