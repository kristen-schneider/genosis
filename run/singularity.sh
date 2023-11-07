#!/usr/bin/env bash

### TODO:
### modify these optinos for your system

#SBATCH -p short
#SBATCH --job-name=example
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32gb
#SBATCH --time=00:20:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/home/name/precision-medicine/example/out/slurm_example.out
#SBATCH --error=/home/name/precision-medicine/example/err/slurm_example.err

set -e pipefail

### TODO: 
### modify these paths

pmed_dir="/home/name/precision-medicine/"
run_script=$pmed_dir"run/example_run.sh"
data_dir=$pmed_dir"example/"
singularity_container="/home/name/pmed.sif"

###
###

bash $run_script


# in one single go -- X
# singularity containers cannot have SLURM installed.
# need to break up the workflow
#singularity run \
#	--bind $data_dir \
#	$singularity_container \
#	bash $run_script
