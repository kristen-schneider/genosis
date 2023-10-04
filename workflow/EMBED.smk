from types import SimpleNamespace
#
config = SimpleNamespace(**config)

shell.prefix("""
. /opt/conda/etc/profile.d/conda.sh
conda activate pmed;
""")

import glob
from os.path import basename

encode_dir=f"{config.out_dir}encodings/"
pos_encodings=glob.glob(encode_dir + "*.pos")
pos_encodings=list(map(basename, pos_encodings))
pos_encodings=[".".join(p.split('.')[:-1]) for p in pos_encodings]
assert len(pos_encodings) > 0, "no positional encodings.."


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
