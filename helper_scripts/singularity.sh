#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=singularity-container
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu
#SBATCH --output=/Users/krsc0813/singularity.out
#SBATCH --error=/Users/krsc0813/singularity.err

pmed_dir="/Users/krsc0813/precision-medicine/"
pmed_sin="/Users/krsc0813/pmed_singularity"

singularity build --sandbox $pmed_sin $pmed_dir
