#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=iLASH-benchmark
#SBATCH --output=/Users/krsc0813/1KG_data/chr15_20/out/iLASH-benchmark.out
#SBATCH --error=/Users/krsc0813/1KG_data/chr15_20/err/iLASH-benchmark.err
#SBATCH --time=0-23:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

# This script shows an example of how we ran iLASH
# including generating ped/map files with plink

# TODO: fill in file paths appropriately.
vcf="/Users/krsc0813/1KG_data/1kg_30x_phased/1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
reference_map="/Users/krsc0813/1KG_data/chr15_20/chrm_15.map"
out_name="/Users/krsc0813/1KG_data/chr15_20/chrm_15.ilash"
ilash_dir="/Users/krsc0813/iLASH/"
ilash_config="/Users/krsc0813/1KG_data/chr15_20/chrm_15.ilash.config"

# create ped/map files from VCF using plink
#plink \
#	--vcf $vcf \
#	--recode 01 \
#	--output-missing-genotype . \
#	--out $out_name

## create map file from interpolated map for iLASH
# awk 'print{$1, $3, $2, $3}' interpolated.map > iLASH.map

# run ilash
$ilash_dir"./build/ilash" $ilash_config
