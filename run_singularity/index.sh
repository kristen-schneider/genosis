#!/usr/bin/env bash

gess_dir=$1
out_dir=$2
config=$3

## These directories should be correct.
## If you have changed where scripts exist, change these paths
smk_dir=$gess_dir"workflow/"
log=$out_dir"log/index.log"
##
##

# for singularity container
. /opt/conda/etc/profile.d/conda.sh
conda activate snakemake

# go to project directory
cd $gess_dir
ls -lah

# run indexing pipeline
# 1. slice vcf
echo "1. slicing VCF..." > $log
#start_slice=$(date +%s.%3N)
snakemake \
    -s $smk_dir"SLICE.smk" \
    -c 16 \
    -j 5 \
    --configfile=$config
#end_slice=$(date +%s.%3N)
#slice_time=$(echo "scale=3; $end_slice - $start_slice" | bc)
#echo "--SLICE: $slice_time seconds" >> $log

# 2. encode slices
echo "2. encode VCF slices..." >> $log
start_encode=$(date +%s.%3N)
snakemake \
    -s $smk_dir"ENCODE.smk" \
    -c 16 \
    -j 10 \
    --configfile=$config \
    --latency-wait 70
#end_encode=$(date +%s.%3N)
#encode_time=$(echo "scale=3; $end_encode - $start_encode" | bc)
#echo "--ENCODE: $encode_time seconds" >> $log

# 3. embed slices
echo "3. embed slice encodings..." >> $log
start_embed=$(date +%s.%3N)
snakemake \
    -s $smk_dir"EMBED_CPU.smk" \
    -c 16 \
    -j 10 \
    --configfile=$config \
    --latency-wait 70
#end_embed=$(date +%s.%3N)
#embed_time=$(echo "scale=3; $end_embed - $start_embed" | bc)
#echo "--EMBED: $embed_time seconds" >> $log

# 4. index slices
echo "4. index slice embeddings..." >> $log
start_index=$(date +%s.%3N)
snakemake \
    -s $smk_dir"INDEX.smk" \
    -c 16 \
    -j 10 \
    --configfile=$config \
    --latency-wait 70
#end_index=$(date +%s.%3N)
#index_time=$(echo "scale=3; $end_index - $start_index" | bc)
#echo "--INDEX: $index_time seconds" >> $log
