#!/bin/env bash

set exuo pipefail

sudo apt-get update
sudo apt install -y awscli tmux vim git

# DONE Test pytorch environment
# DONE Nvidia driver on ML-in-a-box is good enough

# CONDA/MAMBA -----------------------------------------------------------------
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda
eval "$($HOME/miniconda/bin/conda shell.bash hook)"
conda init
conda config --add channels conda-forge bioconda
conda install -y -c conda-forge mamba


# TMUX/VIM --------------------------------------------------------------------
echo "source-file ~/.tmux.d/.tmux.conf" > ~/.tmux.conf
git clone https://github.com/mchowdh200/.tmux.d.git ~/.tmux.d

git clone https://github.com/mchowdh200/.vim.git ~/.vim
echo "export EDITOR=vim" >> ~/.profile

# SNAKEMAKE -------------------------------------------------------------------
mamba env create -f ../envs/snakemake.yml

# TODO 




# TODO post install - wandb login
# TODO post install - copy over data
# TODO post install - edit config or add a commented out version of config that uses cloud machine specific info


