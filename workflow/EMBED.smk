from types import SimpleNamespace
#
config = SimpleNamespace(**config)

#shell.prefix("""
#. /opt/conda/etc/profile.d/conda.sh
#conda activate pmed;
#""")

shell.prefix("""
set -e pipefail
source  ~/.bashrc
conda activate pmed;
""")


import glob
from os.path import basename

rule all:
    input:
        f"{config.out_dir}embeddings/all.embeddings.txt"

# 1.0 run model 
rule model:
    input:
        f"{config.out_dir}zeros.out"
    output:
        f"{config.out_dir}embeddings/all.embeddings.txt"
    message:
        "Running model to create embedding vectors..."
    shell:
        "conda activate torch-gpu;"
        "test ! -d {config.out_dir}embeddings/ && mkdir {config.out_dir}embeddings/;"
        "python {config.root_dir}pytorch/encode_samples.py" \
        " --encoder {config.root_dir}last.ckpt" \
        " --output {config.out_dir}embeddings/all.embeddings.txt" \
        " --gpu" \
        " --files {config.out_dir}encodings/*.pos" \
        " --batch-size {config.batch_size}" \
        " --num-workers {config.n_workers};"
