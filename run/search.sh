#!/usr/bin/env bash

set -e pipefail

pmed_dir=$1
out_dir=$2
config=$3

## These directories should be correct.
## If you have changed where scripts exist, change these paths
smk_dir=$pmed_dir"workflow/"
log=$out_dir"log/search.log"
##
##

# for singularity container
#. /opt/conda/etc/profile.d/conda.sh
# load conda and activate snakemake env for run
#module load anaconda
#.~/miniconda3/bin/conda
#conda_dir="/home/sdp/miniconda3/envs/"
#source ~/miniconda3/etc/profile.d/mamba.sh 
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/miniconda3/condabin/conda
#conda activate snakemake

# activate conda / mamba
source  ~/.bashrc
# go to project directory
cd $pmed_dir

# run searching pipeline
# 5. search slices
echo "5. searching index slices..." >> $log
start_search=$(date +%s.%3N)
snakemake \
    -s $smk_dir"SEARCH.smk" \
    -c 16 \
    -j 10 \
    --configfile=$config \
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
    --configfile=$config \
end_aggregate=$(date +%s.%3N)
aggregate_time=$(echo "scale=3; $end_aggregate - $start_aggregate" | bc)
echo "--AGGREGATE: $aggregate_time seconds" >> $log
