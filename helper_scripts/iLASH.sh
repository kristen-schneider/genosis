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
vcf="/Users/krsc0813/1KG_data/chr15_20/chrm_15-20.vcf.gz"
reference_map="/Users/krsc0813/1KG_data/chr15_20/chrm_15-20.map"
out_name="/Users/krsc0813/1KG_data/chr15_20/chrm_15-20.ilash"
ilash_config="/Users/krsc0813/1KG_data/chr15_20/ilash.config"

# create ped/map files from VCF using plink
plink \
	--vcf $vcf \
	--recode 01 \
	--output-missing-genotype . \
	--out $out_name

# run ilash
./build/ilash $ilash_config
