#!/bin/bash

# VariantViewR_run_pipeline.sh

# Call with one positional argument which will be a prefix for the current study.
# bash VariantViewR_run_pipeline.sh prefix

# Script runs cluster-related steps for VariantViewR pipeline
# Ensure parent directory contains raw_data/ with R1 and R2 fastq.gz files

# Immediately exit on errors
set -ueo pipefail

# Generate lookup.txt file from R1 files in prefix/raw_data/
echo "Writing lookup file"
find $1/raw_data/ -name "*R1*" -exec basename {} \; > $1_lookup.txt

echo "Running array job"
sbatch --wait ./VariantViewR_array.sbatch $1

echo "Concatenating alignment counts info"
# Concatenate a tab-delim header with the output of all the ./sample_data/alignment_data/*_alignCounts.txt
cat <(echo -e "sample\ttotal_reads\tprimary_alignments") \
        ./$1/sample_data/alignment_data/*alignCounts.txt > ./$1/alignment_counts.txt

echo "Done"
