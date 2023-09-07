#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=conda-env
#SBATCH --output=./out/conda-env.out
#SBATCH --error=./err/conda-env.err
#SBATCH --time=0-2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

conda_environment_dir="/Users/krsc0813/precision-medicine/conda_environments/"

eval "$(conda shell.bash hook)"
conda env create -f $conda_environment_dir"pmed.yaml"
conda env create -f $conda_environment_dir"snakemake.yaml"
