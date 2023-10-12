#!/usr/bin/env bash

set -e pipefail

### 
### 

pmed_dir=$1
out_dir=$2
config=$3

smk_dir=$pmed_dir"workflow/"
log=$out_dir"pipeline.log"

##
##

# go to project directory
cd $pmed_dir

# load conda and activate snakemake env for run
. /opt/conda/etc/profile.d/conda.sh
conda activate snakemake

# 1. slice vcf
echo "1. slicing VCF..." > $log
start_slice=$(date +%s.%3N)
snakemake \
    -s $smk_dir"SLICE.smk" \
    -c 16 \
    -j 5 \
    --configfile=$config \
end_slice=$(date +%s.%3N)
slice_time=$(echo "scale=3; $end_slice - $start_slice" | bc)
echo "--SLICE: $slice_time seconds" >> $log

# 2. encode slices
echo "2. encode VCF slices..." >> $log
start_encode=$(date +%s.%3N)
snakemake \
    -s $smk_dir"ENCODE.smk" \
    -c 16 \
    -j 10 \
    --configfile=$config \
end_encode=$(date +%s.%3N)
encode_time=$(echo "scale=3; $end_encode - $start_encode" | bc)
echo "--ENCODE: $encode_time seconds" >> $log
