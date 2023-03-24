#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed

# Given a VCF and a list of sample IDs, this script
# will generate a new VCF only including the
# sample IDs in the list.

# TODO: fill in file paths appropriately.
all_samples_vcf=""
subset_samples_list=""
subset_samples_vcf=""

bcftools view -S $subset_samples_list $all_samples_vcf > $subset_samples_vcf
bgzip -c $subset_samples_vcf
tabix -p vcf $subset_samples_vcf".gz"
