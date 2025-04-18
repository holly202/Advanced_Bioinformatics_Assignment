#!/bin/bash

#6 Using an alternative tool

#Pre-Alignment QC
#Perform quality assessment and trimming
#Move the raw fastq files to untrimmed_fastq directory and rename the files
mv NGS0001.R1.fastq.qz ~/7BBG2016/data/untrimmed_fastq/NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz ~/7BBG2016/data/untrimmed_fastq/NGS0001.R2.fastq.gz
#Quality  assessment
cd ~/7BBG2016/data/untrimmed_fastq
fastqc -t 4 *.fastq.gz
#Store  the results
mkdir  ~/7BBG2016/results/fastqc_untrimmed_reads
mv *fastqc* ~/7BBG2016/results/fastqc_untrimmed_reads/
#Unzip  the other output
cd  ~/7BBG2016/results/fastqc_untrimmed_reads/
for zip in *.zip; do unzip $zip; done
#Cat all fastqc summary.txt files into one fastqc_summaries.txt and move this to  ~/7BBG2016/logs/
cat */summary.txt > ~/7BBG2016/logs/fastqc_summaries.txt

#Trimming 
cd ~/7BBG2016/data/untrimmed_fastq
trimmomatic PE -threads 4  -phred33 ~/7BBG2016/data/untrimmed_fastq/NGS0001.R1.fastq.gz ~/7BBG2016/data/untrimmed_fastq/NGS0001.R2.fastq.gz -baseout ~/7BBG2016/data/trimmed_fastq/NGS0001_trimmed_R ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50
#Perform basic quality assessment of paired trimmed sequencing data
fastqc -t 4 ~/7BBG2016/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/7BBG2016/data/trimmed_fastq/NGS0001_trimmed_R_2P
#Store the results
mkdir ~/7BBG2016/results/fastqc_trimmed_reads
mv ~/7BBG2016/data/trimmed_fastq/*fastqc* ~/7BBG2016/results/fastqc_trimmed_reads/

#Alignment 
genome hg19
#Create a folder for the reference and its index files
mkdir -p ~/7BBG2016/data/reference
mv ~/7BBG2016/data/hg19.fa.gz ~/7BBG2016/data/reference/

#Generate the index files
sudo apt install bowtie2
cd ~/7BBG2016/data/reference
bowtie2-build -f hg19.fa hg19_index
#Replace the aligner with Bowtie2
bowtie2 -x hg19_index -1 ~/7BBG2016/data/trimmed_fastq/NGS0001_trimmed_R_1P -2 ~/7BBG2016/data/trimmed_fastq/NGS0001_trimmed_R_2P -S ~/7BBG2016/data/aligned_data/NGS0001.sam

#Perform duplicate marking
cd ~/7BBG2016/data/aligned_data
#Convert  the sam file into bam format
samtools view -h -b NGS0001.sam > NGS0001.bam
#Sort the bam file
samtools sort NGS0001.bam > NGS0001_sorted.bam
#Generate an index
samtools index NGS0001_sorted.bam file
#Use picard to mark duplicated reads
#The new bam file is stored with duplicate reads marked, the txt file summarises the number of duplicate reads found.
picard MarkDuplicates I=NGS0001_sorted.bam O=NGS0001_sorted_marked.bam M=marked_dup_metrics.txt
samtools index NGS0001_sorted_marked.bam
#Quality Filter the duplicate marked BAM file
#Set Minimum MAPQ quality score : 20
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam
samtools index NGS0001_sorted_filtered.bam
#Generate standard alignment statistics
#Fagstats 
samtools flagstat NGS0001_sorted_marked.bam > NGS0001_sorted_flagstat.txt
samtools flagstat NGS0001_sorted_filtered.bam > NGS0001_filtered_flagstat.txt
#Idxstats
samtools idxstats NGS0001_sorted_filtered.bam > NGS0001_filtered_idxstats.txt
#Insert size
picard CollectInsertSizeMetrics I=NGS0001_sorted_filtered.bam O=insert_size_metrics.txt H=insert_size_histogram.pdf M=0.05
#Depth of coverage
#report  the per-base of coverage for each feature in the annotation.bed
bedtools coverage -a NGS0001_sorted_filtered.bam -b ~/7BBG2016/data/annotation.bed > count_overlaps.bed
#Sort data on column 1,2 and 5
sort -k1,1 -k2,2 -k5,5 count_overlaps.bed > count_overlaps_sorted.bed
#Compute the average depth of coverage
datamash -H mean 6 median 6 min 6 max 6 sstdev 6 < count_overlaps_sorted.bed > depth_coverage.txt

#Variant Calling
#Call Variants using Freebayes
#Convert to text format the reference
zcat ~/7BBG2016/data/reference/hg19.fa.gz > ~/7BBG2016/data/reference/hg19.fa
samtools faidx ~/7BBG2016/data/reference/hg19.fa
#Call variants with Freebayes
freebayes --bam ~/7BBG2016/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/7BBG2016/data/reference/hg19.fa --vcf ~/7BBG2016/results/NGS0001.vcf
#Restrict the VCF to the regions in the bed file provided
bedtools intersect -header -wa -a ~/7BBG2016/results/NGS0001.vcf -b ~/7BBG2016/data/annotation.bed > ~/7BBG2016/results/NGS0001_restricted.vcf
#Compress the resulting variant file
bgzip  ~/7BBG2016/results/NGS0001_restricted.vcf
#Index the VCF
tabix -p vcf ~/7BBG2016/results/NGS0001_restricted.vcf.gz
#Quality Filter Variants
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/7BBG2016//results/NGS0001_restricted.vcf.gz > ~/7BBG2016//results/NGS0001_restricted_filtered.vcf
bgzip  ~/7BBG2016/results/NGS0001_restricted_filtered.vcf
tabix -p vcf ~/7BBG2016/results/NGS0001_restricted_filtered.vcf.gz

#Variant Annotation and Prioritization
#Annotate variants
#Using snpEFF
unzip snpEff_v4_3t_core.zip
java -jar ~/snpEff/snpEff.jar download -v hg19
java -jar ~/snpEff/snpEff.jar hg19 ~/7BBG2016/results/NGS0001_restricted_filtered.vcf > ~/7BBG2016/results/NGS0001_restricted_filtered.ann.vcf
#Using ANNOVAR
bgzip  ~/7BBG2016/results/NGS0001_restricted_filtered.vcf
tabix -p vcf ~/7BBG2016/results/NGS0001_restricted_filtered.vcf.gz
~/annovar/convert2annovar.pl -format vcf4 ~/7BBG2016/results/NGS0001_restricted_filtered.vcf.gz > ~/7BBG2016/results/NGS0001_restricted_filtered.avinput
#Quality Filter Variants
#Use SnpSift filter  
java -jar ~/snpEff/SnpSift.jar filter "(ANN[*].EFFECT has 'exon') && (dbSNP_ID == '')" ~/7BBG2016/results/NGS0001_restricted_filtered.ann.vcf > ~/7BBG2016/results/NGS0001_restricted_filtered_prioritization.vcf
