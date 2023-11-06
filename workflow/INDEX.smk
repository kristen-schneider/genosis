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

# get a list of all positional encodings
ENCODE_DIR=f"{config.out_dir}encodings/"
POS_ENCODINGS=glob.glob(ENCODE_DIR + "*.pos")
POS_ENCODINGS=list(map(basename, POS_ENCODINGS))
POS_ENCODINGS=[pos_enc.replace('.pos', '') for pos_enc in POS_ENCODINGS]
assert len(POS_ENCODINGS) > 0, "no positional encodings.."

rule all:
    input:
        expand(f"{config.out_dir}embeddings/{{segment}}.emb", segment=POS_ENCODINGS),
        svs_config=expand(f"{config.out_dir}svs_index/{{segment}}.config/", segment=POS_ENCODINGS)
        #svs_data=expand(f"{config.out_dir}svs_index/{{segment}}.data/", segment=POS_ENCODINGS),
        #svs_graph=expand(f"{config.out_dir}svs_index/{{segment}}.graph/", segment=POS_ENCODINGS)
        #f"{config.out_dir}svs_index/idx.done"

# 1.0 split all_embeddings.txt into segment embeddings for input to indexing steps
rule split_embeddings:
    input:
        f"{config.out_dir}embeddings/all.embeddings.txt"
    output:
        embeddings=expand(f"{config.out_dir}embeddings/{{segment}}.emb", segment=POS_ENCODINGS)
    message:
        "Splitting full embedding file into segments..."
    shell:
        "conda activate pmed;"
        "python {config.root_dir}python/scripts/split_embeddings.py" \
        " --emb_dir {config.out_dir}embeddings/" \
        " --all_emb {config.out_dir}embeddings/all.embeddings.txt;"

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
        embeddings=f"{config.out_dir}embeddings/{{segment}}.emb"
    output:
        svs_config=directory(f"{config.out_dir}svs_index/{{segment}}.config/"),
        #svs_data=directory(f"{config.out_dir}svs_index/{{segment}}.data/"),
        #svs_graph=directory(f"{config.out_dir}svs_index/{{segment}}.graph/")
    message:
        "SVS-building indexes..."
    shell:
        "conda activate svs;"
        "test ! -d {config.out_dir}svs_index/ && mkdir {config.out_dir}svs_index/;" \
        "python {config.root_dir}python/scripts/svs/build_svs_index.py" \
        " --emb_file {input.embeddings}" \
        " --idx_dir {config.out_dir}svs_index/" \
        " --db_samples {config.out_dir}database_hap_IDs.txt;"
