#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed
#conda list

# running snakemake file
cd /path/to/git-repo/precision-medicine/
git submodule init
git submodule update
echo "...running snakemake"
#snakemake --unlock
snakemake -c
