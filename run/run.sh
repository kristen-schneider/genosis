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
#SBATCH --output=example/log/example.out
#SBATCH --error=example/err/example.err

set -e pipefail

### TODO: 
### modify these paths
# the gess repo
gess_dir="/Users/krsc0813/precision-medicine/"
# the local directory for output and config files
local_dir="/Users/krsc0813/precision-medicine/example/"

### scripts
index_script=$gess_dir"run/index.sh"
search_script=$gess_dir"run/search.sh"
### config file
config_file=$local_dir"example.yml"
###
###

# make log directory
test ! -d $local_dir"log/" && mkdir $local_dir"log/"

# move to project dir
cd $gess_dir

# run indexing step
echo "Running indexing step."
bash $index_script $gess_dir $local_dir $config_file
# run searching step
wait

echo "Running searching step."
bash $search_script $gess_dir $local_dir $config_file
#wait

## run evaluation step
#echo "Running evaluation step."
#bash $evaluate_script $gess_dir $out_dir $config_file
