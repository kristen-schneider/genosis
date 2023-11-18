#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=iLASH-prep
#SBATCH --output=/scratch/Users/krsc0813/chr1_22/out/ilash-prep.out
#SBATCH --error=/scratch/Users/krsc0813/chr1_22/err/ilash-prep.err
#SBATCH --time=0-23:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

# This script shows an example of how we generate files for running iLASH

# TODO: fill in file paths appropriately.
out_dir="/scratch/Users/krsc0813/chr1_22/iLASH/"
vcf_dir="/Users/krsc0813/1KG_data/1kg_30x_phased/"
iLASH_dir="/Users/krsc0813/iLASH/"

for CHROM in {1..22}
do
    echo $CHROM
    config_file=$out_dir"/configs/chrm"$CHROM.yml   

    # running iLASH
    echo "...running iLASH..."
    $ilash_dir"./build/ilash" $config_file
done
