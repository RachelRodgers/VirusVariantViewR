#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=rachel.rodgers

eval "$(conda shell.bash hook)"
conda activate snakemake_7.22.0

snakemake --profile slurm_vvr
