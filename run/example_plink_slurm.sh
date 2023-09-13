#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=plink
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32gb
#SBATCH --time=00:20:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=name@email.com
#SBATCH --output=/home/name/precision-medicine/example/out/plink_example.out
#SBATCH --error=/home/name/precision-medicine/example/err/plink_example.err

set -e pipefail

pmed_dir="/home/name/precision-medicine/"
run_dir=$pmed_dir"run/"
data_dir=$pmed_dir"example/"
singularity_container="/home/name/plink.sif"

singularity run \
	--bind $data_dir \
	$singularity_container \
	bash $run_dir"example_plink_run.sh"
