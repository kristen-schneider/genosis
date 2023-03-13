#!/bin/env bash
#SBATCH --partition=aa100
#SBATCH --gres=gpu:1
#SBATCH --job-name=gt_similarity_search
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=100GB
#SBATCH --time=24:00:00 # hrs:min:sec
#SBATCH --output=log/pairings_pipeline-%j.log
#SBATCH --error=log/pairings_pipeline-%j.err

snakemake -s training_pipeline.smk -j 4 -c 32 --use-conda --conda-frontend mamba
