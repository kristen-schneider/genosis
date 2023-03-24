#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed

# TODO: fill in directories appropriately.
thousand_genome_vcf_dir=""
thousand_genome_map_dir=""
example_dir=""
example_vcf_dir=""

# MERGE VCF
# this script subsets the first 500,000 basepair positions (not variants!)
# from chromosomes 1 through 4, and stores them as separate, small (subset)  vcf files
# then merges those smaller vcf files into 1 larger vcf file
# 1. subset vcfs into smaller vcfs
for i in {1..4}
do
	full_vcf=$thousand_genome_vcf_dir"1kGP_high_coverage_Illumina.chr${i}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
	small_vcf=$example_vcf_dir"chr${i}.small.vcf"
	echo $small_vcf
	
	bcftools view -h $full_vcf > $small_vcf
	tabix $full_vcf chr${i}:0-2000000 >> $small_vcf
	bgzip $small_vcf
	tabix -p vcf $small_vcf".gz"
done
# 2. merge small vcfs
merged_vcf_file=$example_dir"example_merge.vcf"
bcftools concat $example_vcf_dir"chr1.small.vcf.gz" $example_vcf_dir"chr2.small.vcf.gz" $example_vcf_dir"chr3.small.vcf.gz" $example_vcf_dir"chr4.small.vcf.gz" > $merged_vcf_file
bgzip $merged_vcf_file
tabix -p vcf $merged_vcf_file".gz"

# MERGE MAP FILES
# this script concatenates map files from chromosomes 1-5 into one larger map file.
# map files are created with Beagle for 1kg GRCh38.
# download at: https://bochet.gcc.biostat.watshington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
merged_map_file=$example_dir"example_merge.map"
cat $thousand_genome_map_dir"plink.chr1.GRCh38.map" > $merged_map_file
cat $thousand_genome_map_dir"plink.chr2.GRCh38.map" >> $merged_map_file
cat $thousand_genome_map_dir"plink.chr3.GRCh38.map" >> $merged_map_file
cat $thousand_genome_map_dir"plink.chr4.GRCh38.map" >> $merged_map_file

echo "Done"
