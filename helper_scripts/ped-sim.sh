#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed

# This script shows an example of how we generated
# simulated data with ped-sim.

# TODO: fill in file paths appropriately.
data_dir=""
ped_sim_dir=""

$ped_sim_dir"ped-sim" \
	-d $data_dir"example.def" \
	-i $data_dir"example.vcf.gz" \
	-m $data_dir"interpolated.map" \
   	-o $data_dir"out-name" \
	--pois \
	--founder_ids \
	--keep_phase \
	--fam \
	--bp \
	--mrca
	#--retain_extra -1 \
	#--intf $data_dir"chr8.interfere.tsv" \
