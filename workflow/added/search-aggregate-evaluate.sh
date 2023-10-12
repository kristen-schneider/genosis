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

# 5. search slices
echo "5. creating knn from slice indexes..." > $log
start_search=$(date +%s.%3N)
snakemake \
    -s $smk_dir"SEARCH.smk" \
    -c 16 \
    -j 5 \
    --configfile=$config \
end_search=$(date +%s.%3N)
search_time=$(echo "scale=3; $end_search - $start_search" | bc)
echo "--SEARCH: $search_time seconds" >> $log

# 6. aggregate slices
echo "6. aggregating samples from slice KNN..." > $log
start_aggregate=$(date +%s.%3N)
snakemake \
    -s $smk_dir"AGGREGATE.smk" \
    -c 16 \
    -j 5 \
    --configfile=$config \
end_aggregate=$(date +%s.%3N)
aggregate_time=$(echo "scale=3; $end_aggregate - $start_aggregate" | bc)
echo "--AGGREGATE: $aggregate_time seconds" >> $log

# 7. evaluate outcome
echo "7. evaluating and plotting results..." > $log
start_evaluate=$(date +%s.%3N)
snakemake \
    -s $smk_dir"EVALUATE.smk" \
    -c 16 \
    -j 5 \
    --configfile=$config \
end_evaluate=$(date +%s.%3N)
evaluate_time=$(echo "scale=3; $end_evaluate - $start_evaluate" | bc)
echo "--EVALUATE: $evaluate_time seconds" >> $log
