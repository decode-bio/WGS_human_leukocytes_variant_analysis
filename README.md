# Objetive

# Tools
* SRA Toolkit
* Trimmomatic
* BWA
* SAMTOOLS
* GATK
* BCFTOOLS
* SNPEFF

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
bwa mem -t 4 -R "@RG\tID:SRR34149094\tPL:ILLUMINA\tSM:SRR34149094" \
/home/sne_desh/Ref/Homo_sapiens_assembly38.fasta \
/home/sne_desh/WGS/P1/SRR34149094_1.fastq  /home/sne_desh/WGS/P1/SRR34149094_2.fastq > aligned_output.sam
```
<img width="1912" height="912" alt="image" src="https://github.com/user-attachments/assets/f3e3f439-825e-4acc-8f4e-87dab4824832" />

Lets see the file contents using samtools before proceeding to further steps
```
samtools view aligned_output.sam | less
```
<img width="1907" height="989" alt="image" src="https://github.com/user-attachments/assets/eaee39a4-bfe9-4f05-b277-f77bbafa76dd" />

Lets see the statistics of our output
```
samtools flagstat aligned_output.sam
```
<img width="1144" height="437" alt="image" src="https://github.com/user-attachments/assets/41d51bf7-649a-43f2-836e-2812cd9083b0" />


Now we get our aligned read in .sam format. Lets continue with the next step to do sorting and duplicate removal.

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
samtools -b -S aligned_output.sam > aligned_unsorted.bam
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
samtools markdup -r aligned_position.fixmate.bam aligned_final.bam
```

Indexing the final sorted bam file
```
samtools index aligned_final.bam
```
Output
<img width="1888" height="124" alt="image" src="https://github.com/user-attachments/assets/034e7305-6f81-4a27-add7-556fed6050c4" />

```
samtools view aligned_final.bam | less
```
<img width="1906" height="953" alt="image" src="https://github.com/user-attachments/assets/4ce89ce2-472e-4aa2-a7da-5084ddec3e8c" />


Now I am done with the sorting and removal of duplicate and indexing of my bam file. lets go to the next lets check the statistics of our aligned reads.

* Statistics of the Aligned reads
```
samtools flagstat aligned_final.bam
```
<img width="1117" height="449" alt="image" src="https://github.com/user-attachments/assets/22588c43-8ff5-457c-ba8f-0e0026dc5bb0" />


Overall mapping is 99.73%, which indicates that our sample matches the reference & library preparedness was clean.

# Step 6: Variant calling using GATK
Now before going Forward wth variant calling, I checked for
* Sample name
* Header
```
samtools view aligned_final.bam | less
```
<img width="1528" height="811" alt="image" src="https://github.com/user-attachments/assets/8b10f7c0-1c8c-4e82-995a-15a0f6f2265d" />

```
samtools view -H aligned_final.bam | grep '^@RG'
```
<img width="1131" height="50" alt="image" src="https://github.com/user-attachments/assets/1f2f2ef4-30f5-46ce-a820-ae03562fa3c4" />

As my ``` aligned_final.bam ``` contains both sample_name and header. So we can  dive into the variant calling using GATK.

Now lets perform the GATK variant calling 
```
python3 ~/Tools/gatk-4.6.2.0/gatk \
  --java-options "-Xmx8g" HaplotypeCaller  \
  -R $ref \
  -I aligned_final.bam \
  -O Output.vcf \
  -ERC GVCF 
```
<img width="1911" height="1010" alt="image" src="https://github.com/user-attachments/assets/60cd343f-0d33-48da-9e79-22dab745530b" />

<img width="1919" height="1015" alt="image" src="https://github.com/user-attachments/assets/01efea88-42cb-4f2a-b45a-5fd4acf67423" />

Now before separating the variants, lets count the variants in the .vcf file using the bcf tools
```
bcftools +counts Output.vcf
```
<img width="1549" height="171" alt="image" src="https://github.com/user-attachments/assets/6fb5776b-d163-47c0-bf8d-5bef112e944f" />
Now lets separate the snps and Indels into two different files. For separating we can either use GATK's Selectvariants or bcftools. Here I am using the bcfstools to separate the variants.

* SNP

```
bcftools view -v snps Output.vcf -O v -o Output_snp.vcf
```

* INDEL

```
bcftools view -v indels Output.vcf -O v -o Output_indels.vcf
```

Step 7: Variant Annotation
Now lets performe the variant annotation using the snpEff tool.
```
snpEff -v -stats report.html hg38 output.vcf>annotated.vcf
```
The code will generate the three files
* annotated.vcf
* report.html
* report.genes.txt
<img width="1912" height="945" alt="image" src="https://github.com/user-attachments/assets/4bcc5515-dc60-4d81-94d1-09cea4776033" />



