#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=rachel.rodgers

ml snakemake/5.10.0-python-3.6.5

snakemake --profile slurm_vvr
