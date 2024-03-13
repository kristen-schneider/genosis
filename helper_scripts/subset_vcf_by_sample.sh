#!/usr/bin/env bash

### TODO:
### modify these options for your system

#SBATCH -p long
#SBATCH --job-name=subset_samples
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32gb
#SBATCH --time=23:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/Users/krsc0813/gess_1KG_pops/log/subset_samples_EUR_21.log
#SBATCH --error=/Users/krsc0813/gess_1KG_pops/err/subset_samples_EUR_21.err

set -e pipefail

# loading modules
#echo "...loading modules."
#module load anaconda

# creating and activating conda environment
#echo "...activating conda environment."
#conda activate pmed

# Given a VCF and a list of sample IDs, this script
# will generate a new VCF only including the
# sample IDs in the list.

# TODO: fill in file paths appropriately.
all_samples_vcf="/Users/krsc0813/1KG_data/1kg_30x_phased/1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
subset_samples_list="/Users/krsc0813/gess_1KG_pops/EUR_samples.txt"
subset_samples_vcf="/Users/krsc0813/gess_1KG_pops/EUR_chr21.vcf"

bcftools view -S $subset_samples_list $all_samples_vcf > $subset_samples_vcf
bgzip $subset_samples_vcf
tabix -p vcf $subset_samples_vcf".gz"
