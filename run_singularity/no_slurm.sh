#!/usr/bin/env bash

set -e pipefail

pmed_dir="/home/sdp/precision-medicine/"

data_dir=$pmed_dir"example/"
#data_dir="/Users/krsc0813/chr10/"
#data_dir="/Users/krsc0813/chr15_20/"
#data_dir="/Users/krsc0813/AFR_pedigree/"

log=$data_dir"pipeline.log"

# go to project directory and update
cd $pmed_dir
git submodule init
git submodule update

eval "$(conda shell.bash hook)"
#conda info --envs

# run pipeline
# 1. slice vcf
echo "1. slicing VCF..." > $log
start_slice=$(date +%s.%3N)
#conda activate snakemake
snakemake \
    -s "./run_singularity/"SLICE.smk \
    -c 16 \
    -j 5 \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete
end_slice=$(date +%s.%3N)
#slice_time=$(echo "scale=3; $end_slice - $start_slice" | bc)
#echo "--SLICE: $slice_time seconds" >> $log

# 2. encode slices
echo "2. encode VCF slices..." >> $log
start_encode=$(date +%s.%3N)
snakemake \
    -s "./run_singularity/"ENCODE.smk \
    -c 16 \
    -j 10 \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete
end_encode=$(date +%s.%3N)
#encode_time=$(echo "scale=3; $end_encode - $start_encode" | bc)
#echo "--ENCODE: $encode_time seconds" >> $log

# 3. embed slices
echo "3. embed slice encodings..." >> $log
conda activate base
start_embed=$(date +%s.%3N)
echo $start_embed
snakemake \
    -s "./run_singularity/"EMBED.smk \
    -c 16 \
    -j 10 \
    --rerun-incomplete
    #--use-conda \
    #--conda-frontend mamba \
    #--cluster-config embed_config.yaml \
    #--cluster "sbatch -J {cluster.job-name} \\
    #                  -t {cluster.time} \\
    #                  -N {cluster.nodes} \\
    #                  -n {cluster.ntasks} \\
    #                  -p {cluster.partition} \\
    #                  --nodelist={cluster.nodelist} \\
    #                  --gres={cluster.gpu} \\
    #                  --mem={cluster.mem} \\
    #                  --output {cluster.output} \\
    #                  --error {cluster.error}" \
    #--latency-wait 70 \
end_embed=$(date +%s.%3N)
echo $end_embed
#embed_time=$(echo "scale=3; $end_embed - $start_embed" | bc)
#echo "--EMBED: $embed_time seconds" >> $log

# 4. index slices
echo "4. index slice embeddings..." >> $log
start_index=$(date +%s.%3N)
snakemake \
    -s "./run_singularity/"INDEX.smk \
    -c 16 \
    -j 10 \
    --use-conda \
    --conda-frontend mamba \
    --rerun-incomplete
end_index=$(date +%s.%3N)
#index_time=$(echo "scale=3; $end_index - $start_index" | bc)
#echo "--INDEX: $index_time seconds" >> $log
#
## 5. search slices
#echo "5. searching index slices..." >> $log
#start_search=$(date +%s.%3N)
#snakemake \
#    -s SEARCH.smk \
#    -c 16 \
#    -j 10 \
#    --use-conda \
#    --conda-frontend mamba \
#    --rerun-incomplete
#end_search=$(date +%s.%3N)
#search_time=$(echo "scale=3; $end_search - $start_search" | bc)
#echo "--SEARCH: $search_time seconds" >> $log
#
#
## 6. aggregate slices
#echo "6. aggregating results slices..." >> $log
#start_aggregate=$(date +%s.%3N)
#snakemake \
#    -s AGGREGATE.smk \
#    -c 16 \
#    -j 10 \
#    --use-conda \
#    --conda-frontend mamba \
#    --rerun-incomplete
#end_aggregate=$(date +%s.%3N)
#aggregate_time=$(echo "scale=3; $end_aggregate - $start_aggregate" | bc)
#echo "--AGGREGATE: $aggregate_time seconds" >> $log
