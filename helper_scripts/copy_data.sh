#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=copy-data
#SBATCH --output=./out/copy-data.out
#SBATCH --error=./err/copy-data.err
#SBATCH --time=0-23:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

origin_dir="/scratch/Users/krsc0813/chr1_22/svs_results/"
new_dir="/scratch/Users/krsc0813/copy_chr1_22/"

mkdir $new_dir
cp -r $origin_dir $new_dir
