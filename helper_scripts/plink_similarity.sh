#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed

# This script shows an example of how we ran
# plink --genome and 
# plink2 --make-king-table

# TODO: fill in file paths appropriately.
vcf=""
data_dir=""

plink --vcf $vcf \
       	--genome \
	--out $data_dir"pink-genome"

plink2 --vcf /path/to/vcf.vcf.gz --make-king-table`
