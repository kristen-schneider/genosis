#!/bin/env bash
#SBATCH --partition=amilan
#SBATCH --job-name=gt_similarity_search
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=64G
#SBATCH --time=2:00:00 # hrs:min:sec
#SBATCH --output=/dev/null
#SBATCH --error=log/pairings_pipeline_%j.err

snakemake -s exploration_pipeline.smk -j 4 --use-conda --conda-frontend mamba
