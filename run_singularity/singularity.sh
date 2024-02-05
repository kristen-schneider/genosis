#!/usr/bin/env bash

### TODO:
### modify these optinos for your system

#SBATCH -p short
#SBATCH --job-name=singularity
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32gb
#SBATCH --time=00:20:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/Users/krsc0813/precision-medicine/example/log/singularity_run.out
#SBATCH --error=/Users/krsc0813/precision-medicine/example/err/singularity_run.err

set -e pipefail

### TODO: 
### modify these paths

pmed_dir="/Users/krsc0813/precision-medicine/"
run_script=$pmed_dir"run_singularity/run.sh"
data_dir=$pmed_dir"example/"
singularity_container="/Users/krsc0813/pmed_singularity/gess.sif"

###
###

#bash $run_script


# in one single go -- X
# singularity containers cannot have SLURM installed.
# need to break up the workflow
singularity run \
	--bind $data_dir \
	$singularity_container \
	bash $run_script
