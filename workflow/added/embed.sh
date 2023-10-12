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

# 3. embed slices
echo "3. creating embeddings from slice encodings..." > $log
start_embed=$(date +%s.%3N)
snakemake \
    -s $smk_dir"EMBED.smk" \
    -c 16 \
    -j 5 \
    --configfile=$config \
end_embed=$(date +%s.%3N)
embed_time=$(echo "scale=3; $end_embed - $start_embed" | bc)
echo "--EMBED: $embed_time seconds" >> $log
