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
#SBATCH --output=/Users/user/example/out/example.out
#SBATCH --error=/Users/user/example/err/example.err

set -e pipefail

### TODO: 
### modify path to run script
run_script="/Users/user/example/example_merge.sh"
###
###

bash $run_script
