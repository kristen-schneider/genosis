#!/usr/bin/env bash

gess_dir=$1
out_dir=$2
config=$3

## These directories should be correct.
## If you have changed where scripts exist, change these paths
smk_dir=$gess_dir"workflow/"
log=$out_dir"log/search.log"
##
##

# for singularity container
. /opt/conda/etc/profile.d/conda.sh
conda activate snakemake

# go to project directory
cd $gess_dir

# run searching pipeline
# 5. search slices
echo "5. searching index slices..." >> $log
start_search=$(date +%s.%3N)
snakemake \
    -s $smk_dir"SEARCH.smk" \
    -c 16 \
    -j 10 \
    --configfile=$config \
    --latency-wait 70
end_search=$(date +%s.%3N)
search_time=$(echo "scale=3; $end_search - $start_search" | bc)
echo "--SEARCH: $search_time seconds" >> $log


# 6. aggregate slices
echo "6. aggregating results slices..." >> $log
start_aggregate=$(date +%s.%3N)
snakemake \
    -s $smk_dir"AGGREGATE.smk" \
    -c 16 \
    -j 10 \
    --configfile=$config
end_aggregate=$(date +%s.%3N)
aggregate_time=$(echo "scale=3; $end_aggregate - $start_aggregate" | bc)
echo "--AGGREGATE: $aggregate_time seconds" >> $log
