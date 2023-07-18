#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=aggregate
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64gb
#SBATCH --time=23:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu
#SBATCH --output=/Users/krsc0813/precision-medicine/slurm_scripts/out/aggregate.out
#SBATCH --error=/Users/krsc0813/precision-medicine/slurm_scripts/err/aggregate.err

set -e pipefail

pmed_dir="/Users/krsc0813/precision-medicine/"

#data_dir="/Users/krsc0813/chr10/"
data_dir="/Users/krsc0813/chr10_12/"
#data_dir=$pmed_dir"example/"
#data_dir="/Users/krsc0813/AFR_pedigree/"

log=$data_dir"pipeline.log"

# go to project directory and update
cd $pmed_dir
git submodule init
git submodule update


# 6. aggregate slices
echo "6. aggregating results slices..." >> $log
start_aggregate=$(date +%s.%3N)
snakemake \
    -s AGGREGATE.smk \
    -c 16 \
    -j 10 \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete
end_aggregate=$(date +%s.%3N)
aggregate_time=$(echo "scale=3; $end_aggregate - $start_aggregate" | bc)
echo "--AGGREGATE: $aggregate_time seconds" >> $log
