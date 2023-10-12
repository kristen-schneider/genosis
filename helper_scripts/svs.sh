#!/usr/bin/env bash

#SBATCH -p short
#SBATCH --job-name=svs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64gb
#SBATCH --time=1:00:00
#SBATCH --mail-type=NONE
#SBATCH --mail-user=krsc0813@colorado.edu
#SBATCH --output=/Users/krsc0813/svs.out
#SBATCH --error=/Users/krsc0813/svs.err

svs_wheel='/Users/krsc0813/pysvs-0.0.1-cp311-cp311-manylinux_2_17_x86_64.manylinux2014_x86_64.whl'


pip install $svs_wheel --force-reinstall
