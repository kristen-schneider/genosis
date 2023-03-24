#!/bin/bash

# loading modules
echo "...loading modules."
module load anaconda

# creating and activating conda environment
echo "...activating conda environment."
conda activate pmed

# TODO: fill in file paths appropriately.
python_file="./python/scripts/ibd/phased_ibd.py"
data_dir=""
vcf_file=""
map_file=$data_dir"interpolated.map"

python $python_file \
	--vcf $vcf_file \
	--map $map_file \
	> $data_dir"example.phasedibd.out"
