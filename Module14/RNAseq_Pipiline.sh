#!/bin/bash

SECONDS=0

# change working directory
cd <path>


# STEP 1: Run fastqc
fastqc demo.fastq 

# run trimmomatic to trim reads with poor quality
java -jar ~/Desktop/demo/tools/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 demo.fastq demo_trimmed.fastq TRAILING:10 -phred33

fastqc demo_trimmed.fastq 



# STEP 2: Run HISAT2
# get the genome indices
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz

# run alignment
hisat2 -q --rna-strandness R -x <path>/grch38/genome -U demo_trimmed.fastq | samtools sort -o demo_trimmed.bam
samtools view demo_trimmed.bam



# STEP 3: Run featureCounts - Quantification
# get gtf
# wget https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.gtf.gz

featureCounts -S 2 -a <path>/Homo_sapiens.GRCh38.106.gtf -o demo_featurecounts.txt demo_trimmed.bam
echo "featureCounts finished running!"

cat demo_featurecounts.txt.summary | less
cat demo_featurecounts.txt | less
cat demo_featurecounts.txt | cut -f1,7 | less




