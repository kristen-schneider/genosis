#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=plot
#SBATCH --output=./out/plot.out
#SBATCH --error=./err/plot.err
#SBATCH --time=0-23:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

plot_script="/Users/krsc0813/biobagg_analysis/plotting/evaluate_ancestry.py"
pop_file="/Users/krsc0813/1KG_data/samples.ancestry"
knn_file="/scratch/Users/krsc0813/chr1_22/TOP_HITS.txt"
png_dir="/scratch/Users/krsc0813/chr1_22/violin_plots/"

python $plot_script \
    --pop $pop_file\
    --knn $knn_file\
    --png $png_dir
