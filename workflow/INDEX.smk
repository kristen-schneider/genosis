from types import SimpleNamespace
#
config = SimpleNamespace(**config)

#shell.prefix("""
#. /opt/conda/etc/profile.d/conda.sh
#conda activate pmed;
#""")

shell.prefix("""
source  ~/.bashrc
conda activate pmed;
""")

import glob
from os.path import basename

EMBED_DIR=f"{config.out_dir}embeddings/"
EMBEDDINGS=glob.glob(EMBED_DIR + "*.emb")
EMBEDDINGS=list(map(basename, EMBEDDINGS))
EMBEDDINGS=[".".join(p.split('.')[:-1]) for p in EMBEDDINGS]
assert len(EMBEDDINGS) > 0, "no embeddings.."

rule all:
    input:
        expand(f"{config.out_dir}embeddings/{{segment}}.emb", segment=EMBEDDINGS),
        f"{config.out_dir}svs_index/idx.done"
	#f"{config.out_dir}log/faiss_build.log",

# 1.0 split all_embeddings.txt into segment embeddings for input to indexing steps
rule split_embeddings:
    input:
        f"{config.out_dir}embeddings/all.embeddings.txt"
    output:
        expand(f"{config.out_dir}embeddings/{{segment}}.emb", segment=EMBEDDINGS),
    message:
        "Splitting full embedding file into segments..."
    shell:
        "conda activate pmed;"
        "python {config.root_dir}python/scripts/split_embeddings.py" \
        " --emb_dir {config.out_dir}embeddings/" \
        " --all_emb {config.out_dir}embeddings/all.embeddings.txt;" \

## 5.1 build FAISS indices
#rule faiss_build:
#    input:
#        emb_segments=f"{config.out_dir}embeddings/{{segment}}.emb"
#    output:
#        idx_segments=f"{config.out_dir}faiss_index/{{segment}}.idx"
#    message:
#        "FAISS-building indexes..."
#    conda:
#        f"{config.conda_faiss}"
#    shell:
#        "echo 4. ---CREATING FAISS INDEX---;" \
#        "test ! -d {config.out_dir}faiss_index/ && mkdir {config.out_dir}faiss_index/;" \
#        "python {config.root_dir}python/scripts/faiss/build_faiss_index.py" \
#        " --emb {input.emb_segments}" \
#        " --idx_dir {config.out_dir}faiss_index/" \
#        " --db_samples {config.database_IDs}" \
#        " --pedigree {config.pedigree}"
#

# 5.2 build SVS indices
rule svs_build:
    input:
        expand(f"{config.out_dir}embeddings/{{segment}}.emb", segment=EMBEDDINGS)
    output:
        idx_done=f"{config.out_dir}svs_index/idx.done"
    message:
        "SVS-building indexes..."
    shell:
        "conda activate base;"
        "test ! -d {config.out_dir}svs_index/ && mkdir {config.out_dir}svs_index/;" \
        "python {config.root_dir}python/scripts/svs/build_svs_index.py" \
        " --emb_dir {config.out_dir}embeddings/" \
        " --idx_dir {config.out_dir}svs_index/" \
        " --db_samples {config.out_dir}database_hap_IDs.txt" \
        " --emb_ext emb;" \
        "touch {output.idx_done};"
