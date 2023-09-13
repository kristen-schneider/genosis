#!/usr/bin/env bash

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

pmed_dir="/home/name/precision-medicine/"
run_dir=$pmed_dir"run/"
data_dir=$pmed_dir"example/"
singularity_container="/home/name/pmed.sif"

singularity run \
	--bind $data_dir \
	$singularity_container \
	bash $run_dir"example_run.sh"
