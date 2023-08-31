#!/usr/bin/env bash

#SBATCH --partition=short
#SBATCH --job-name=download-vcf
#SBATCH --output=./out/download-vcf.out
#SBATCH --error=./err/download-vcf.err
#SBATCH --time=0-2:00:00
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu

vcf_dir="/Users/krsc0813/1kg_30x_phased/"
vcf_file="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
vcf_file_idx="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"

cd $vcf_dir
wget $vcf_file
wget $vcf_file_idx

