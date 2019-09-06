# VirusVariantViewR
------
VariantViewR consists of a data analysis pipeline and R Shiny application originally developed for the analysis of viral variant information (specifically, murine norovirus CR6 variants). The application continues to be modified and developed.

The full pipeline consists of several steps prior to interaction with the application.

The recommended directory structure is a parent directory containing the following subdirectories:
* parent_directory/VirusVariantViewR/ (scripts found in this repository)
* parent_directory/Mod_CR6_ORFs/
* subdirectory for each data set to be run through the VariantViewR pipeline scripts (step 1 below) with the raw fastq data residing within a raw_data/ subdirectory (ie: parent_dir/study_directory/raw_samples).  Additional subdirectories will be added to the study_directory/ when running the pipeline scripts.

## Data Preparation 

1. **VariantViewR_run_pipeline.sh** (which calls **VariantViewR_array.sbatch**) processes raw sequencing data and generates VCF files.  This pipeline consists of the following steps:
	1. adapter and quality trimming
	2. deduplication to remove exact (PCR) duplicates
	3. alignment to CR6 reference genome with bowtie2
	4. variant calling with bcftools mpileup and bcftools calling
	5. generation of additional information of interest including calculating average coverage and number of primary alignments per sample.
	
These scripts will generate a specific directory structure that should be preserved when running the mutation annotation R script and the Shiny app:

study_directory/
	alignment_files/
	QC_reads/
	raw_data/ (user-generated)
	sample_data/
	variants/
	alignment_counts.txt
	

2. **Mod_CR6_ORF_Information.R** reads files found in the Mod_CR6_ORFs/ directory and writes the R data object "Mod_CR6_ORF_Data.RData."  Please ensure the Mod_CR6_ORFs/ directory is one level above the script, or modify the path to these files as needed. The generated R data object is then read into the script **Annotation_Mutations.R**.
	
3. **Annotate_Mutations.R** annotates each mutation contained the VCF files found in the study_directory/variants/ directory.  These annotated files are required for the shiny application.
	
## Data Display/Interaction with VariantViewR App

4. Application requires **global.R** and **app.R**.  Application is launched from **app.R**.

