# VirusVariantViewR
------
VariantViewR consists of a data analysis pipeline and R Shiny application originally developed for the analysis of viral variant information (specifically, murine norovirus CR6 variants). The application continues to be modified and developed.

The full pipeline consists of several steps prior to interaction with the application.

## Data Preparation 

1. Processing of raw data to VCF files with the scripts **VariantViewR_run_pipeline.sh** (which calls **VariantViewR_array.sbatch**).  This pipeline consists of the following steps:
	1. adapter and quality trimming
	2. deduplication to remove exact (PCR) duplicates
	3. alignment to CR6 reference genome with bowtie2
	4. variant calling with bcftools mpileup and bcftools calling
	5. generation of additional information of interest including calculating average coverage and number of primary alignments per sample.
	
2. Annotation of the mutations contained in each VCF file with the script **Annotate_Mutations.R**.  These annotated files are required for the shiny application.
	
3. Generation of "Mod_CR6_ORF_Information.RData" required for the shiny application with the script **Mod_CR6_ORF_Information.R**. 

## Data Display/Interaction with VariantViewR App

4. Application requires **global.R** and **app.R**.  Application is launched from **app.R**.

