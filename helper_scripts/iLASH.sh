#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed

# This script shows an example of how we ran iLASH
# including generating ped/map files with plink

# TODO: fill in file paths appropriately.
vcf=""
reference_map=""
out_name=""
ilash_config=""

# create ped/map files from VCF using plink
plink \
	--vcf $vcf \
	--recode 01 \
	--output-missing-genotype . \
	--out $out_name

# run ilash
./build/ilash $ilash_config
