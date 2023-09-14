#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=example
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32gb
#SBATCH --time=00:20:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/home/kristens/kristen-pmed/decode/out/slurm_decode.out
#SBATCH --error=/home/kristens/kristen-pmed/decode/err/slurm_decode.err

set -e pipefail

### TODO: 
### modify these paths

pmed_dir="/home/kristens/kristen-pmed/precision-medicine/"
run_script="/home/kristens/kristen-pmed/decode/decode_run.sh"
data_dir="/home/kristens/kristen-pmed/decode/"
singularity_container="/home/kristens/kristen-pmed/pmed_singularity/pmed.sif"

###
###

singularity run \
	--bind $data_dir \
	$singularity_container \
	bash $run_script
