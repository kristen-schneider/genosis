from types import SimpleNamespace
#
config = SimpleNamespace(**config)

shell.prefix("""
. /opt/conda/etc/profile.d/conda.sh
conda activate pmed;
""")

import glob
from os.path import basename

ENCODE_DIR=f"{config.out_dir}encodings/"
POS_ENCODINGS=glob.glob(ENCODE_DIR + "*.pos")
POS_ENCODINGS=list(map(basename, POS_ENCODINGS))
POS_ENCODINGS=[".".join(p.split('.')[:-1]) for p in POS_ENCODINGS]
assert len(POS_ENCODINGS) > 0, "no positional encodings.."


rule all:
    input:
        f"{config.out_dir}embeddings/model.out",

# 1.0 run model 
rule model:
    input:
        f"{config.out_dir}zeros.out"
    output:
        model_out=f"{config.out_dir}embeddings/model.out"
    resources:
        slurm_partition="fijigpu-04"
    message:
        "Running model to create embedding vectors..."
    shell:
        "conda activate torch-cpu;"
        "test ! -d {config.out_dir}embeddings/ && mkdir {config.out_dir}embeddings/;" \
        "python {config.root_dir}pytorch/encode_samples.py" \
            " --encoder {config.root_dir}last.ckpt" \
            " --output {config.out_dir}embeddings/all.embeddings.txt" \
        " --gpu" \
        " --files {config.out_dir}encodings/*.pos" \
        " --batch-size {config.batch_size}" \
        " --num-workers {config.n_workers};" \
	"touch {output.model_out};"
