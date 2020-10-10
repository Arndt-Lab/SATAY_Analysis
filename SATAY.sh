#!/bin/bash
# see http://core.sam.pitt.edu/node/5678
#SBATCH -N 1
#SBATCH -t 0-08:00:00 # 8hrs
#SBATCH -c 8 # Request 4 cpus-per-task. #8
#SBATCH --mem=128g 
#SBATCH --job-name=SATAY


#set aliases
# file_to_analyze_prefix=[YOUR FASTQ FILE NAME PREFIX]
# 
# dir=/bgfs/karndt/[YOUR SATAY ANALYSIS DIRECTORY]
# fastq_dir=/bgfs/karndt/[YOUR SATAY FASTQ DIRECTORY]
# satay_index=/bgfs/karndt/satay_index/satay_index

#What File do you want to analyze?
file_to_analyze_prefix=WT_DpnII

dir=/bgfs/karndt/SATAY_TESTING_MITCH
fastq_dir=/bgfs/karndt/SATAY_BP/fastq
satay_index=/bgfs/karndt/satay_index/satay_index
chromosome_sizes=/bgfs/karndt/satay_index

#Do you want stats on all FASTQ and BAM FIles???
#YES means you do and NO means you dont
get_stats=YES
# get_stats=NO

#Which Enzyme did you digest with??
enzyme=DpnII
# enzyme=NlaIII


##LOAD MODULES
module load deeptools/3.3.0
module load prinseq/0.20.4
module load cutadapt/2.10
module load kentutils/3.0.2 #ucsc utilities
module load bowtie2/2.4.1
module load gcc/8.2.0
module load bedtools/2.29.0
module load samtools/1.9


#SETUP OUTPUT DIRECTORY STRUCTURE
#set working directory path
cd ${dir}

#For folder removal when we want to rerun the whole analysis from scratch
# rm -rf trimmed_fastq
# rm -rf trimmed_fastq_stats
# rm -rf alignment_stats
# rm -rf bam
# rm -rf bam_stats
# rm -rf tag_bed
# rm -rf tag_bigwig
# rm -rf tag_counts

#re/create folders for data
mkdir unzipped_fastq
mkdir unzipped_fastq_stats
mkdir trimmed_fastq
mkdir trimmed_fastq_stats
mkdir alignment_stats
mkdir bam
mkdir bam_stats
mkdir full_bigwig
mkdir tag_bed
mkdir tag_bigwig
mkdir read_counts_all_genes
mkdir read_counts_essential_genes


#TRIM READS WITH CUTADAPT
#.fastq.gz is a compressed/zipped form of the .fastq file. .fastq is the raw output of the Illumina sequencer.
#need to unzip the files first
gunzip -c ${fastq_dir}/${file_to_analyze_prefix}.fastq.gz > ${dir}/unzipped_fastq/${file_to_analyze_prefix}.fastq

#Calculate FASTQ Stats
if [ ${get_stats} = "YES" ]
then
  prinseq-lite.pl -stats_all -fastq ${dir}/unzipped_fastq/${file_to_analyze_prefix}.fastq > ${dir}/unzipped_fastq_stats/${file_to_analyze_prefix}_stats.txt
fi



#This removes any transposon sequence from the 3' end of reads so that the ONLY sequence mapping is insert sequence.
#This also throws out reads that are smaller than 20 bp to better uniquely map reads.

if [ ${enzyme} = "DpnII" ]
then
	#THIS COMMAND IS SPECIFIC TO LIBRARIES BUILD FROM DpnII DIGESTED FRAGMENTS
	cutadapt -a GATCGTATCGGTTTTCGATTACCGTATTTATCCCGTTCGTTTTCGT -m 20 -o \
		${dir}/trimmed_fastq/${file_to_analyze_prefix}_trimmed.fastq \
		${dir}/unzipped_fastq/${file_to_analyze_prefix}.fastq
fi

if [ ${enzyme} = "NlaIII" ]
then
	#THIS COMMAND IS SPECIFIC TO LIBRARIES BUILD FROM NlaIII DIGESTED FRAGMENTS
	cutadapt -a CATGTGTGATTTTACCGAACAAAAATACCGGTTCCCGTCCGATTTCGACTTTAACCCGACCGGATCGTATCGGTTTTCGATTACCGTATTTATCCCGTTCGTTTTCGT -m 20 -o \
		${dir}/trimmed_fastq/${file_to_analyze_prefix}_trimmed.fastq \
		${dir}/unzipped_fastq/${file_to_analyze_prefix}.fastq
fi


#Calculate FASTQ Stats
if [ ${get_stats} = "YES" ]
then
prinseq-lite.pl -stats_all -fastq ${dir}/trimmed_fastq/${file_to_analyze_prefix}_trimmed.fastq > ${dir}/trimmed_fastq_stats/${file_to_analyze_prefix}_trimmed_stats.txt
fi


#ALIGN READS WITH BOWTIE2
bowtie2 --no-unal \
	-q \
	-k 2 \
	-p 4 \
	-x ${satay_index} \
	-U ${dir}/trimmed_fastq/${file_to_analyze_prefix}_trimmed.fastq \
	2> ${dir}/alignment_stats/${file_to_analyze_prefix}_bowtie2_stats.txt \
	| samtools view - -bS -q 30 -F 0x0400 -@ 4 > ${dir}/bam/${file_to_analyze_prefix}_unsorted.bam

#SORT BAM
samtools sort -l 9 -O BAM -@ 4 -T ${dir}/sc_bam/${dir}/bam/${file_to_analyze_prefix}_sorted \
	-o ${dir}/bam/${file_to_analyze_prefix}_sorted.bam \
	${dir}/bam/${file_to_analyze_prefix}_unsorted.bam
	
	
#INDEX BAM FOR VIEWING IN IGV
samtools index -@ 4 ${dir}/bam/${file_to_analyze_prefix}_sorted.bam


#Calculate FASTQ Stats
if [ ${get_stats} = "YES" ]
then
samtools flagstat -@ 4 ${dir}/bam/${file_to_analyze_prefix}_sorted.bam >  ${dir}/bam_stats/${file_to_analyze_prefix}_flagstat.txt
samtools stats -@ 4 ${dir}/bam/${file_to_analyze_prefix}_sorted.bam | grep ^SN | cut -f 2- >  ${dir}/bam_stats/${file_to_analyze_prefix}_stats.txt
samtools idxstats -@4 ${dir}/bam/${file_to_analyze_prefix}_sorted.bam | cut -f 1,3 >  ${dir}/bam_stats/${file_to_analyze_prefix}_chr_and_pdna_alnstats.txt
fi


#GENERATE BROWSER TRACKS

#full reads
bamCoverage --bam ${dir}/bam/${file_to_analyze_prefix}_sorted.bam \
	--outFileName ${dir}/full_bigwig/${file_to_analyze_prefix}_full_reads.bw \
	--outFileFormat bigwig \
	--binSize 1 \
	--numberOfProcessors 4

#convert to bedgraph file with 5' tag counts
bedtools genomecov -bg -5 -ibam ${dir}/bam/${file_to_analyze_prefix}_sorted.bam > ${dir}/tag_bed/${file_to_analyze_prefix}.bedgraph

#convert to bigwig
bedGraphToBigWig -blockSize=1 ${dir}/tag_bed/${file_to_analyze_prefix}.bedgraph ${chromosome_sizes}/chr_size.txt ${dir}/tag_bigwig/${file_to_analyze_prefix}.bw


#COUNT READS OVER GENES
#all
bedtools coverage \
	-a /bgfs/karndt/Annotations/BED_files/mRNA_sorted_fixed.bed \
	-b ${dir}/bam/${file_to_analyze_prefix}_sorted.bam \
	-counts \
	> ${dir}/read_counts_all_genes/${file_to_analyze_prefix}_read_counts.tsv
#essential
bedtools coverage \
	-a /bgfs/karndt/Annotations/BED_files/essential_genes_fix.bed \
	-b ${dir}/bam/${file_to_analyze_prefix}_sorted.bam \
	-counts \
	> ${dir}/read_counts_essential_genes/${file_to_analyze_prefix}_read_counts_essential.tsv





exit