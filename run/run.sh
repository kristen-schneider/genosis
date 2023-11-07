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

pmed_dir="./"
index_script="./run/index.sh"
search_script="./run/search.sh"
data_dir="./example/"

###
###

cd $pmed_dir
# run indexing step
echo "Running indexing step."
bash $index_script
echo "Running searching step."
bash $search_script
