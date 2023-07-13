#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=clean
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu
#SBATCH --output=/Users/krsc0813/precision-medicine/slurm_scripts/out/clean.out
#SBATCH --error=/Users/krsc0813/precision-medicine/slurm_scripts/err/clean.err

data_dir='/Users/krsc0813/chr10/'
#data_dir='/Users/krsc0813/chr10_12/'
#data_dir='/Users/krsc0813/AFR_pedigree/'
#data_dir='/Users/krsc0813/precision-medicine/example/'

rm /Users/krsc0813/precision-medicine/slurm_scripts/err/chr10.err
rm /Users/krsc0813/precision-medicine/slurm_scripts/out/chr10.out

# go to data directory 
cd $data_dir

rm *.txt
rm segment_boundary.map
rm interpolated.map

rm -r log/
rm -r benchmark/
rm -r vcf_segments/
rm -r encodings/
rm -r embeddings/
rm -r faiss_index/
rm -r faiss_results/
rm -r query_results/
