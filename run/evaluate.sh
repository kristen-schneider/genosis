#!/usr/bin/env bash

set -e pipefail

### TODO: 
### modify these paths

pmed_dir="./"
out_dir="./example/"
config=$out_dir"example.yml"

###
###

## These directories should be correct.
## If you have changed where scripts exist, change these paths
smk_dir=$pmed_dir"workflow/"
log=$out_dir"log/evaluate.log"
##
##

# go to project directory
cd $pmed_dir

# for singularity container
#. /opt/conda/etc/profile.d/conda.sh
# load conda and activate snakemake env for run
#module load anaconda
#.~/miniconda3/bin/conda
#conda_dir="/home/sdp/miniconda3/envs/"
#source ~/miniconda3/etc/profile.d/mamba.sh 
echo $SHELL
source  ~/.bashrc
#export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/miniconda3/condabin/conda
#conda activate snakemake

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
