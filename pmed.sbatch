#!/bin/bash
# submits snakemake jobs as individual slurm jobs
# -----------------------------------------------
#SBATCH --partition ami100
##SBATCH --gres=gpu
#SBATCH --job-name=pmed-snakemake
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=0-0:30:00
#SBATCH --qos=normal
#SBATCH --output=/scratch/alpine/krsc0813/slurm_files/out/pmed.out
#SBATCH --error=/scratch/alpine/krsc0813/slurm_files/err/pmed.err
#SBATCH --mail-type=None
#SBATCH --mail-user=krsc0813@colorado.edu
# -----------------------------------------------



# activate conda env (slurm)
echo "Preparing conda env."
echo "...loading anaconda module."
module load anaconda
echo "...activating conda environment."
conda activate pmed
# activate conda env (bash)
#conda_source='/home/sdp/miniconda3/etc/profile.d/conda.sh'
#source $conda_source
#conda activate pmed


# environment variables
project_dir='/scratch/alpine/krsc0813/precision-medicine/'
#project_dir='/home/sdp/precision-medicine/'
sm_config=$project_dir'example/config_cluster.yaml'

# update packages
echo "Configuring libraries."
cd /scratch/alpine/krsc0813/precision-medicine/
git submodule init
git submodule update


# run snakemake
echo "Running snakemake."
snakemake -s Snakefile \
        --use-conda --conda-frontend mamba \
        --cluster-config $sm_config \
        --latency-wait 5 \
        --cluster "sbatch -J {cluster.job-name} \\
                          -t {cluster.time} \\
                          -N {cluster.nodes} \\
                          -n {cluster.ntasks} \\
                          --error {cluster.error}" \
        --jobs 4 --cores 1
