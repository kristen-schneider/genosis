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

# 4. index slices
echo "4. creating indexes from slice embeddings..." > $log
start_index=$(date +%s.%3N)
snakemake \
    -s $smk_dir"INDEX.smk" \
    -c 16 \
    -j 5 \
    --configfile=$config \
end_index=$(date +%s.%3N)
index_time=$(echo "scale=3; $end_index - $start_index" | bc)
echo "--INDEX: $index_time seconds" >> $log
