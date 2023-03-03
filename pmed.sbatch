#!/bin/bash
# submits snakemake jobs as individual slurm jobs
# -----------------------------------------------
#SBATCH -p short
#SBATCH --job-name=pmed-snakemake
#SBATCH --mail-type=None (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=krsc0813@colorado.edu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=1gb
#SBATCH --time=168:00:00 # Time limit hrs:min:sec
#SBATCH --output=/Users/much8161/Repositories/seq-similarity-search/log/bwa-snakemake.out
#SBATCH --error=/Users/much8161/Repositories/seq-similarity-search/log/bwa-snakemake.err
# -----------------------------------------------



# activate conda env
conda_source='/home/sdp/miniconda3/etc/profile.d/conda.sh'
source $conda_source
conda activate pmed


# environment variables
project_dir='/home/sdp/precision-medicine/'
sm_config=$project_dir'example/config_snakemake.yaml'


# run snakemake
snakemake -c \
	--use-conda --conda-frontend mamba \
	--cluster-config $sm_config \
	--latency-wait 5 \
	--cluster "sbatch -J {cluster.job-name} \\
			  -t {cluster.tme} \\
			  -N {cluster.nodes} \\
			  -n {cluster.ntasks} \\
			  --mem={cluster.mem} \\
			  --mem={cluster.mem} \\
			  --error {cluster.error}" \
	--jobs 4 --cores 64
