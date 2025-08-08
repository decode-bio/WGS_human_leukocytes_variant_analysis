# Objetive

# Tools
* SRA Toolkit
* Trimmomatic
* BWA
* SAMTOOLS
* GATK

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

* On comparing both the reports, it is clearly seen that the sequence base lenght has been reduced, which will increase the false positive alignment. So I will proceed further with the untrimmed fastq files for the downstream analysis. Lets align sequence with the reference genome using BWA-MEM.

# Step4: Alignment using BWA
Now after QC, we proceed with our original reads and run alignment using BWA. For this step reference genome assembly is already been downloaded and indexed.
```
bwa mem  /home/sne_desh/Ref/Homo_sapiens_assembly38.fasta  /home/sne_desh/WGS/P1/SRR34149094_1.fastq  /home/sne_desh/WGS/P1/SRR34149094_2.fastq > aligned.sam
```
<img width="1587" height="959" alt="image" src="https://github.com/user-attachments/assets/1c093ecb-38da-4132-b376-6517109c1ec2" />
<img width="1896" height="718" alt="image" src="https://github.com/user-attachments/assets/83d45197-b3a9-4d8f-99ef-7ae77df9f275" />


Our Overall alignement took total 2.77 hrs. Now we get our aligned read in .sam format. Lets continue with the next step to do sorting and duplicate removal.

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
samtools -b -S aligned.sam > aligned_unsorted.bam
```
Sorting the unsorted .bam file
```
samtools sort -n -o aligned_sorted.bam aligned_unsorted.bam
```
Fix mate-pair information
```
samtools fixmate -m aligned_sorted.bam aligned_sorted.fixmate.bam
```
Sort by genomic coordinates
```
samtools sort -o aligned_position.fixmate.bam aligned_sorted.fixmate.bam
```
Mark and remove duplicates
```
samtools markdup -r aligned_position.fixmate.bam aligned_sorted.rm.bam
```

Indexing the final sorted bam file
```
samtools index aligned_sorted.rm.bam
```
Output
<img width="1671" height="152" alt="image" src="https://github.com/user-attachments/assets/59cefcb5-8943-458f-a984-6d1d77fb3010" />

Now I am done with the sorting and removal of duplicate and indexing of my bam file. lets go to the next lets check the statistics of our aligned reads.

* Statistics of the Aligned reads
```
samtools flagstat aligned_sorted.rm.bam
```
<img width="1114" height="437" alt="image" src="https://github.com/user-attachments/assets/bb62a48b-20eb-4fcc-bc51-d8bec5e7acc4" />

Overall mapping is 99.73%, which indicates that our sample matches the reference & library preparedness was clean.

# Step 6: Variant calling using GATK
Now before going Forward wth variant calling, I checked for
* Sample name
* Header
As my ``` aligned_sorted.bam ``` neither contains sample_name nor header. So before diving into the variant calling the following steps has been performed.

```
python3 ~/Tools/gatk-4.6.2.0/gatk AddOrReplaceReadGroups \
-I aligned_sorted.rm.bam \
-O output_RG.bam \
    -RGID H0164.2 \
    -RGLB library1 \
    -RGPL illumina \
    -RGPU H0164ALXX140820.2 \
    -RGSM sample1

```
Now have to indexed the output file again
```
samtools index output_RG.bam
```
Now lets perform the GATK variant calling 
```
python3 ~/Tools/gatk-4.6.2.0/gatk \
  --java-options "-Xmx8g" HaplotypeCaller  \
  -R $ref \
  -I output_RG.bam \
  -O Output.vcf \
  -ERC GVCF 
```
<img width="1911" height="469" alt="image" src="https://github.com/user-attachments/assets/a222cf8f-5bf3-4be5-958b-92fe2185347b" />
<img width="1810" height="999" alt="image" src="https://github.com/user-attachments/assets/8a0f6550-84a2-4154-870b-01fc859767b2" />
<img width="1910" height="1012" alt="image" src="https://github.com/user-attachments/assets/84b30525-cd1b-48d5-baf7-f102aba77029" />







