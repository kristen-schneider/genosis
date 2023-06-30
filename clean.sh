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

#data_dir='/Users/krsc0813/chr10/'
data_dir='/Users/krsc0813/precision-medicine/example/'

# go to data directory 
cd $data_dir

rm *.txt
rm *.map

rm -r vcf_segments/
rm -r encodings/
rm -r embeddings/
rm -r faiss_index/
rm -r faiss_results
