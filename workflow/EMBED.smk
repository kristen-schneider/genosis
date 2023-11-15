from types import SimpleNamespace
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

# get a list of all encoding segments
ENC_DIR=f"{config.out_dir}encodings/"
ENC_SEGMENTS=glob.glob(ENC_DIR + "*.gt")
ENC_SEGMENTS=list(map(basename, ENC_SEGMENTS))
ENC_SEGMENTS=[enc_seg.replace('.gt', '') for enc_seg in ENC_SEGMENTS]
assert len(ENC_SEGMENTS) > 0, "no genotype encoding segments.."

rule all:
    input:
        all_embeddings=f"{config.out_dir}embeddings/all.embeddings.txt"

# 1.0 generate embedding vectors
rule model:
    input:
        zeros=expand(f"{config.out_dir}encodings/{{segment}}.done", segment=ENC_SEGMENTS)
    output:
        f"{config.out_dir}embeddings/all.embeddings.txt"
    message:
        "Running model to create embedding vectors..."
    shell:
        """
        conda activate torch-gpu;
        test ! -d {config.out_dir}embeddings/ && mkdir {config.out_dir}embeddings/;
        python {config.root_dir}pytorch/encode_samples.py\
         --encoder {config.root_dir}last.ckpt\
         --output {config.out_dir}embeddings/all.embeddings.txt\
         --gpu\
         --files {config.out_dir}encodings/*.pos\
         --batch-size {config.batch_size}\
         --num-workers {config.n_workers};
        """
