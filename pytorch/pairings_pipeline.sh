#!/bin/env bash
#SBATCH -j gt_similarity_search
#SBATCH -p amilan
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=64G
#SBATCH --time=2:00:00 # hrs:min:sec
#SBATCH --output=log/pairings_pipeline_%j.log
#SBATCH --error=log/pairings_pipeline_%j.err

snakemake -s exploration_pipeline.smk -j 4 --use-conda --conda-frontend mamba
