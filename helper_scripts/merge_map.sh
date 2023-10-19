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


# concat map files from plink maps
merged_file="/Users/krsc0813/1KG_data/chr1_22/chrm_1-22.map"
touch $merged_file

cd /Users/krsc0813/plink_maps/


for CHROM in {1..22}
do
    #echo plink.chr$CHROM."GRCh38.map"
    cat "./plink.chr"$CHROM."GRCh38.map" >> $merged_file
done

#cat plink.chr15.GRCh38.map > $merged_file
#cat plink.chr16.GRCh38.map >> $merged_file
#cat plink.chr17.GRCh38.map >> $merged_file
#cat plink.chr18.GRCh38.map >> $merged_file
#cat plink.chr19.GRCh38.map >> $merged_file
#cat plink.chr20.GRCh38.map >> $merged_file
