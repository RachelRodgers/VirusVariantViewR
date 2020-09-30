#!/bin/bash

#----- calculate_alignment_data.sh -----#
#
# For use with VariantViewR Snakefile
#
# Usage:
# bash calculate_alignemnt_data.sh \
#	sample \
#	sample_deduplicated_R1.fastq \
#	sample_deduplicated_R2.fastq \
#	sample_sorted.bam \
#	sample_alignCounts.txt \
#
# Rachel Rodgers, Sep 2020
#---------------------------------------#

# Count QC'd reads

echo "Counting total reads."

read1_total_reads=$(cat $2 | echo $((`wc -l`/4)))
read2_total_reads=$(cat $3 | echo $((`wc -l`/4)))
sample_total_reads=$((read1_total_reads + read2_total_reads))

# Count the number of primary alignments from the sorted bam file
# Will need to add left mate mapped and right mate mapped, then
# dividing this sum by sample_total_reads will yield percent MNV

echo "Counting primary alignments."

left_mate_mapped=$(samtools view -f 0x40 -F 0x4 $4 | cut -f1 | sort | uniq | wc -l)
right_mate_mapped=$(samtools view -f 0x80 -F 0x4 $4 | cut -f1 | sort | uniq | wc -l) 

# Send this information to a new text file

echo -e $1"\t"${sample_total_reads}"\t"${left_mate_mapped}"\t"${right_mate_mapped} > $5