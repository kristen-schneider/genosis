#!/bin/bash

#SBATCH --partition=short
#SBATCH --job-name=merge-map
#SBATCH --output=./out/merge-map.out
#SBATCH --error=./err/merge-map.err
#SBATCH --time=0-01:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu


# concat vcf files from 1000 genomes
merged_file="/Users/krsc0813/chr15_20/chrm_15-20.map"

cd /Users/krsc0813/plink_maps/

cat plink.chr15.GRCh38.map > $merged_file
cat plink.chr16.GRCh38.map >> $merged_file
cat plink.chr17.GRCh38.map >> $merged_file
cat plink.chr18.GRCh38.map >> $merged_file
cat plink.chr19.GRCh38.map >> $merged_file
cat plink.chr20.GRCh38.map >> $merged_file
