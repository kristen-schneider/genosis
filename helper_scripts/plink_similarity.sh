#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=plink-benchmark
#SBATCH --output=/Users/krsc0813/1KG_data/chrm15/out/plink-benchmark.out
#SBATCH --error=/Users/krsc0813/1KG_data/chrm15/err/plink-benchmark.err
#SBATCH --time=0-23:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

# This script shows an example of how we ran
# plink --genome and 
# plink2 --make-king-table

# TODO: fill in file paths appropriately.
vcf="/Users/krsc0813/1KG_data/1kg_30x_phased/1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
data_dir="/Users/krsc0813/1KG_data/chrm15/"

cd $data_dir
echo "plink --genome"
plink --vcf $vcf \
       	--genome \
	--out $data_dir"/plink-genome/chrm15"

echo "plink2 --kinship"
plink2 --vcf $vcf --make-king-table
