#!/usr/bin/env bash

### TODO:
### modify these options for your system

#SBATCH -p short
#SBATCH --job-name=example
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32gb
#SBATCH --time=23:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=./example/out/example.out
#SBATCH --error=./example/err/example.err

set -e pipefail

### TODO: 
### modify these paths
# scripts
index_script="/Users/krsc0813/precision-medicine/run/index.sh"
search_script="/Users/krsc0813/precision-medicine/run/search.sh"
evaluate_script="/Users/krsc0813/precision-medicine/run/evaluate.sh"
# directories
pmed_dir="/Users/krsc0813/precision-medicine/"
out_dir="/Users/krsc0813/precision-medicine/example/"
# config file
config_file=$out_dir"example.yml"
###
###

# make log directory
test ! -d $out_dir"log/" && mkdir $out_dir"log/"

# move to project dir
cd $pmed_dir

# run indexing step
echo "Running indexing step."
bash $index_script $pmed_dir $out_dir $config_file
# run searching step

echo
echo "Running searching step."
bash $search_script $pmed_dir $out_dir $config_file

echo
## run evaluation step
#echo "Running evaluation step."
#bash $evaluate_script $pmed_dir $out_dir $config_file
