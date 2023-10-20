#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=merge-vcf
#SBATCH --time=0-24:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu
#SBATCH --output=/Users/krsc0813/1KG_data/chr1_22/out/merge-vcf.out
#SBATCH --error=/Users/krsc0813/1KG_data/chr1_22/err/merge-vcf.err

# concat vcf files from 1000 genomes
merged_file="/Users/krsc0813/1KG_data/chr1_22/chrm_1-22.vcf"

cd "/Users/krsc0813/1KG_data/1kg_30x_phased/"

bcftools concat \
        1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr3.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr4.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr6.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr7.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr8.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr11.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr12.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr13.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr14.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr16.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr17.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr18.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
        > $merged_file

bgzip $merged_file
tabix -p vcf $merged_file".gz"
