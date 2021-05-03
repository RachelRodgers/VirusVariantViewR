# VirusVariantViewR
------
Snakemake workflow for VirusVariantViewR pipeline

Please visit the VirusVariantViewR application:

https://baldridge-lab-wustl.shinyapps.io/VirusVariantViewR/

------

VariantViewR consists of a data analysis pipeline (Snakemake worfkow) and R Shiny application originally developed for the analysis of viral variant information (specifically, murine norovirus CR6 variants). The pipeline and application also work with murine norovirus CW3 variants. The application continues to be modified and developed.

The full pipeline consists of several steps prior to interaction with the application.

This repository is for the *workflow*. Please visit [here](https://github.com/RachelRodgers/VirusVariantViewR-RShiny-Application) for the companion R Shiny application.

This Snakemake workflow can be used on different computing systems.  Currently it is written for the HTCF cluster using the SLURM manager and a Snakemake profile. Running on a different machine may require the deletion of certain "module load" or "ml" commands within the rules of the Snakefile, and editing the config file(s) appropriately.

How to use:

1. Clone the repository (for HTCF users, clone to your /scratch/ directory):
```
git clone --recurse-submodules https://github.com/RachelRodgers/VirusVariantViewR.git
```
2. For HTCF or other Slurm users, make a directory to hold the snakemake profile:
```
mkdir -p ~/.config/snakemake/slurm_vvr
```
3. Edit the config.yaml file inside the config directory to your user email:
```
--mail-user=yourEmailAddress
```
4. Copy the cluster submit and profile files to the appropriate locations:
```
cd VirusVariantViewR
cp config/config.yaml ~/.config/snakemake/slurm_vvr
cp slurm-submit/*.py ~/.config/snakemake
```
5. Create a directory to hold your raw sequencing reads and move your data into that directory.
6. Edit the vvr_config.yaml file (under /config/) to adjust the data set name, path to your reads, and other parameters as needed.
7. Submit in one of two ways:
	
	a. With sbatch script (preferred for HTCF/Slurm) users:
	```
	sbatch submit_vvr_snake.sbatch
	```
	b. Interactively (better for troubleshooting):
	```
	# start an interactive session (for HTCF users):
	srun --mem=48G --cpus-per-task=8 -J interactive -p interactive --pty /bin/bash -l
	
	# load snakemake:
	ml snakemake/5.10.0-python-3.6.5
	
	# dry run (prints steps and stops):
	snakemake --profile slurm_vvr -n
	
	# production run (run steps):
	snakemake --profile slurm_vvr
	```
(HTCF/Slurm users): See Slurm output files in logs_slurm/ directory which will generate inside the hecatomb_htcf_snake/ directory.

8. Pipeline output will be in the results directory (contained within a directory matching your DataSetName).
9. Visualize the data locally with the [VirusVariantViewR-RShiny-Application](https://github.com/RachelRodgers/VirusVariantViewR-RShiny-Application)!
