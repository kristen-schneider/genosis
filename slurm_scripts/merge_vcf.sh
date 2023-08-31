#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=merge-vcf
#SBATCH --output=./out/merge-vcf.out
#SBATCH --error=./err/merge-vcf.err
#SBATCH --time=0-24:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

# concat vcf files from 1000 genomes
merged_file="/Users/krsc0813/chr15_20/chrm_15-20.vcf"

cd "/Users/krsc0813/1kg_30x_phased/"

bcftools concat \
        1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr16.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr17.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr18.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        > $merged_file

bgzip $merged_file
tabix -p vcf $merged_file".gz"
