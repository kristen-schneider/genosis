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

for CHROM in {1..22}
do
    vcf_file=$vcf_dir"1kGP_high_coverage_Illumina.chr"$CHROM".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    echo $vcf_file
    config_file=$out_dir"/configs/chrm"$CHROM.yml   

    # create PED file
    echo "...creating PED file..."
    plink \
        --vcf $vcf_file \
        --recode 01 \
        --output-missing-genotype . \
        --out $out_dir"chrm"$CHROM
    
    # create MAP file from interpolated map
    echo "...creating MAP file..."
    grep "chrm"$CHROM $out_dir"iLASH.fullmap" > $out_dir"chrm"$CHROM".map"
done

