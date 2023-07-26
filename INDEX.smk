from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/AFR/AFR_config.yaml"
configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/Users/krsc0813/chr10/config_fiji.yaml"
#configfile: "/Users/krsc0813/chr10_12/config_snakemake.yaml"
#configfile: "/Users/krsc0813/AFR_pedigree/config_AFR.yaml"

config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

import glob
from os.path import basename

emb_dir=f"{config.embeddings_dir}"
emb_segments=glob.glob(emb_dir + "*.emb")
emb_segments=list(map(basename, emb_segments))
emb_segments=[".".join(e.split('.')[:-1]) for e in emb_segments]
assert len(emb_segments) > 0, "no embeddings.."

rule all:
    input:
        expand(f"{config.faiss_index_dir}{{segment}}.idx", segment=emb_segments),
	#f"{config.log_dir}faiss_build.log",

# 5.1 build FAISS indices
rule faiss_build:
    input:
        emb_segments=f"{config.embeddings_dir}{{segment}}.emb"
    output:
        idx_segments=f"{config.faiss_index_dir}{{segment}}.idx"
    message:
        "FAISS-building indexes..."
    conda:
        f"{config.conda_faiss}"
    shell:
        "echo 4. ---CREATING FAISS INDEX---;" \
        "test ! -d {config.faiss_index_dir} && mkdir {config.faiss_index_dir};" \
        "python {config.python_dir}faiss/build_faiss_index.py" \
        " --emb {input.emb_segments}" \
        " --idx_dir {config.faiss_index_dir}" \
        " --db_samples {config.database_IDs}" \
        " --pedigree {config.pedigree}"
