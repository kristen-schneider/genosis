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
#SBATCH --output=/Users/krsc0813/precision-medicine/example/log/example.out
#SBATCH --error=/Users/krsc0813/precision-medicine/example/err/example.err

set -e pipefail

### TODO: 
### modify these paths
# scripts
index_script="/Users/krsc0813/precision-medicine/run/index.sh"
search_script="/Users/krsc0813/precision-medicine/run/search.sh"
# directories
gess_dir="/Users/krsc0813/precision-medicine/"
out_dir="/Users/krsc0813/precision-medicine/example/"
# config file
config_file=$out_dir"example.yml"
###
###

# make log directory
test ! -d $out_dir"log/" && mkdir $out_dir"log/"

# move to project dir
cd $gess_dir

# run indexing step
echo "Running indexing step."
bash $index_script $gess_dir $out_dir $config_file
# run searching step
wait

echo "Running searching step."
bash $search_script $gess_dir $out_dir $config_file
#wait

## run evaluation step
#echo "Running evaluation step."
#bash $evaluate_script $gess_dir $out_dir $config_file
