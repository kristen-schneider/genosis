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
#SBATCH --output=/Users/krsc0813/precision-medicine/example/out/example.out
#SBATCH --error=/Users/krsc0813/precision-medicine/example/err/example.err

set -e pipefail

### TODO: 
### modify path to run script
root_dir="/Users/krsc0813/precision-medicine/example/"
run_script=$root_dir"example_merge.sh"
###
###

cd $root_dir
bash $run_script
