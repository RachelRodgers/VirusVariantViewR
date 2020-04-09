# VirusVariantViewR
------

Please visit the VirusVariantViewR application:

[https://baldridge-lab-wustl.shinyapps.io/]

------

VariantViewR consists of a data analysis pipeline and R Shiny application originally developed for the analysis of viral variant information (specifically, murine norovirus CR6 variants). The pipeline and application also work with murine norovirus CW3 variants.  The application continues to be modified and developed.

The full pipeline consists of several steps prior to interaction with the application.

The recommended directory structure is a parent directory containing the following subdirectories:
* parent_directory/VirusVariantViewR/ (scripts found in this repository)
* parent_directory/Mod_CR6_ORFs/
* parent_directory/Mod_CW3_ORFs/
* subdirectory for each data set to be run through the VariantViewR pipeline scripts (step 1 below) with the raw fastq data residing within a raw_data/ subdirectory (ie: your_parent_dir/your_study_directory/raw_samples).  Additional subdirectories will be added to your_study_directory/ when running the pipeline scripts.
* parent_directory/VirusVariantViewR_datasets.txt:  this is a single-column plain text file containing the name of the subdirectory for each data set that has been run through the VariantViewR pipeline scripts, and for which you want to visualize in the application.
* parent_directory/<name_of_dataset>_ metadata.txt: a tab-delimited text file containing your sample names in column 1 (named "Sample") and an arbitary number of named columns with additional data of intereste for each sample.

## Data Preparation 

1. **VariantViewR_run_pipeline.sh** (which calls **VariantViewR_array.sbatch**) processes raw sequencing data and generates VCF files.  This pipeline consists of the following steps:
	1. adapter and quality trimming
	2. deduplication to remove exact (PCR) duplicates
	3. alignment to CR6 reference genome with bowtie2 (requires bt2-indexed reference genome; CR6 reference genome contained in the file /Mod_CR6_ORFs/"Mod_CR6_FullGenome.txt")
	4. variant calling with LoFreq or bcftools mpileup/call
	5. generation of additional information of interest including calculating average coverage and number of primary alignments per sample.
	
These scripts will generate a specific directory structure that should be preserved when running the mutation annotation R script and the Shiny app:

* your_study_directory/alignment_files/
* your_study_directory/QC_reads/
* your_study_directory/raw_data/ (user-generated)
* your_study_directory/sample_data/
* your_study_directory/variants/
* your_study_directory/alignment_counts.txt

2. **Mod_CR6_ORF_Information.R** (if comparing against CR6 genome) reads files found in the Mod_CR6_ORFs/ directory and writes the R data object "Mod_CR6_ORF_Data.RData."  Please ensure the Mod_CR6_ORFs/ directory is one level above the script, or modify the path to these files as needed. The generated R data object is then read into the script **Annotation_Mutations.R**.

3. **Mod_CW3_ORF_Information.R** (if comparing against CW3 genome), same as step 2.
	
3. **Annotate_Mutations.R** annotates each mutation contained the VCF files found in the study_directory/variants/ directory.  These annotated files are required for the shiny application.
	
## Data Display/Interaction with VariantViewR App

4. Application requires **global.R** and **app.R**.  Application is launched from **app.R**.

