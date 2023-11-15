#!/usr/bin/env bash

set -e pipefail

pmed_dir=$1
out_dir=$2
config=$3

## These directories should be correct.
## If you have changed where scripts exist, change these paths
smk_dir=$pmed_dir"workflow/"
log=$out_dir"log/evaluate.log"
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

# run evaluation pipeline
# 7. evaluate KNN
echo "7. evaluating sample results..." >> $log
start_evaluate=$(date +%s.%3N)
snakemake \
    -s $smk_dir"EVALUATE.smk" \
    -c 16 \
    -j 10 \
    --configfile=$config \
end_evaluate=$(date +%s.%3N)
evaluate_time=$(echo "scale=3; $end_evaluate - $start_evaluate" | bc)
echo "--EVALUATE: $evaluate_time seconds" >> $log
