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
merged_file="/Users/krsc0813/chr10_12/chrm_10-12.map"

cd /Users/krsc0813/plink_maps/

cat plink.chr10.GRCh38.map > $merged_file
cat plink.chr11.GRCh38.map >> $merged_file
cat plink.chr12.GRCh38.map >> $merged_file
