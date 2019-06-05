#!/bin/bash

# VariantViewR_array.sbatch

#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --array=1-19%10
#SBATCH -J Variant_ViewR_array
#SBATCH --mail-user=rachel.rodgers
#SBATCH --mail-type=ALL

# Immediately exit on errors
set -ueo pipefail

module load bbtools
module load samtools
module load bowtie2
module load bcftools
module load fastqc
module load bedtools

# Make directories if needed
mkdir -p ./$1/QC_reads/fastQC
mkdir -p ./$1/alignment_files
mkdir -p ./$1/variants
mkdir -p ./$1/sample_data/alignment_data
mkdir -p ./$1/sample_data/genome_coverage

# Fetch read1 from lookup.txt, generate read2
# lookup.txt lives in the same directory as script
read1=$(sed -n ${SLURM_ARRAY_TASK_ID}p $1_lookup.txt)
read2=$(echo ${read1} | sed s/R1/R2/g)

# From the read1 file, get the short sample identifier, which will be of the form: Baldridge_##
# Baldridge_## will be the first two fields of the file name if delimited on underscores
#sample=$(echo ${read1} | cut -d'_' -f 1-2)
sample=$(echo ${read1} | cut -d'_' -f 1-3)

echo "Working on ${sample}."

#----- QC -----#

# adapter trimming
echo "Trimming adapters."
bbduk.sh \
	in1=./$1/raw_data/${read1} \
	in2=./$1/raw_data/${read2} \
	out1=./$1/QC_reads/${sample}_adaptTrim_R1.fastq \
	out2=./$1/QC_reads/${sample}_adaptTrim_R2.fastq \
	ref=~/contaminants_11.20.18/Baldridge_nextera.fa \
	ktrim=r k=23 mink=11 hdist=1 tpe tbo &> ./$1/QC_reads/${sample}_adapterTrim.log

# quality trimming
echo "Quality trimming."
bbduk.sh \
	in1=./$1/QC_reads/${sample}_adaptTrim_R1.fastq \
	in2=./$1/QC_reads/${sample}_adaptTrim_R2.fastq \
	out1=./$1/QC_reads/${sample}_qualTrim_R1.fastq \
	out2=./$1/QC_reads/${sample}_qualTrim_R2.fastq \
	qtrim=r trimq=25 minlength=75 &> ./$1/QC_reads/${sample}_qualityTrim.log

#----- Count QC'd Reads -----#

echo "Counting total reads."
read1_total_reads=$(cat ./$1/QC_reads/${sample}_qualTrim_R1.fastq | echo $((`wc -l`/4)))
read2_total_reads=$(cat ./$1/QC_reads/${sample}_qualTrim_R2.fastq | echo $((`wc -l`/4)))
sample_total_reads=$((read1_total_reads + read2_total_reads))

#----- FastQC QC'd Reads -----#

fastqc ./$1/QC_reads/${sample}_qualTrim_R1.fastq \
	./$1/QC_reads/${sample}_qualTrim_R2.fastq \
	-o ./$1/QC_reads/fastQC

#----- Align to CR6 with Bowtie2 -----#

echo "Aligning to CR6."

# Note - not writing out any unaligned sequences.
# samtools flagstat on sorted bam output will count only mapped reads!

bowtie2 \
	--time \
	--no-unal \
	-x ~/mnv_genomes_11.20.18/Mod_CR6 \
	-1 ./$1/QC_reads/${sample}_qualTrim_R1.fastq \
	-2 ./$1/QC_reads/${sample}_qualTrim_R2.fastq \
	-S ./$1/alignment_files/${sample}_alignment.sam &> ./$1/alignment_files/${sample}_alignment.log

#----- Convert, Sort, Index Alignment Files -----#

echo "Sorting and indexing."

samtools view -S -b ./$1/alignment_files/${sample}_alignment.sam > ./$1/alignment_files/${sample}_alignment.bam
samtools sort ./$1/alignment_files/${sample}_alignment.bam -o ./$1/alignment_files/${sample}_sorted.bam
samtools index ./$1/alignment_files/${sample}_sorted.bam ./$1/alignment_files/${sample}_sorted.bam.bai

echo "Making bedGraph."
# Generate a bedGraph file for the sorted bam file
bedtools genomecov -bga \
        -ibam ./$1/alignment_files/${sample}_sorted.bam > \
        ./$1/alignment_files/${sample}_sorted.bedGraph

echo "Counting primary alignments."

# Count the number of primary alignments from the sorted bam file
# Dividing this number by sample_total_reads will yield percent MNV
sample_primary_alignments=$(samtools view -F 0x4 ./$1/alignment_files/${sample}_sorted.bam | \
	cut -f 1 | sort | uniq | wc -l)

# Send this information to a new text file
echo -e ${sample}"\t"${sample_total_reads}"\t"${sample_primary_alignments} > \
	./$1/sample_data/alignment_data/${sample}_alignCounts.txt
#----- Call Variants -----#

# calling variants relative to Modified_CR6.fas

echo "Calling variants."

# --annotate AD adds (unfiltered) allelic depth - ie: number of reads supporting each allele - to output VCF file

bcftools mpileup \
	-f ~/mnv_genomes_11.20.18/Modified_CR6.fas \
	./$1/alignment_files/${sample}_sorted.bam \
	--annotate AD | \
	bcftools call -m --variants-only > ./$1/variants/${sample}_variants.vcf

# Calculate the average depth of coverage
bedtools genomecov -ibam \
	./$1/alignment_files/${sample}_sorted.bam > \
	./$1/sample_data/genome_coverage/${sample}_coverage.txt

echo "Done with ${sample}."
