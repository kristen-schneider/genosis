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

ENCODE_DIR=f"{config.out_dir}encodings/"
POS_ENCODINGS=glob.glob(ENCODE_DIR + "*.pos")
POS_ENCODINGS=list(map(basename, POS_ENCODINGS))
POS_ENCODINGS=[".".join(p.split('.')[:-1]) for p in POS_ENCODINGS]
assert len(POS_ENCODINGS) > 0, "no positional encodings.."

rule all:
    input:
        embeddings=f"{config.out_dir}embeddings/{{segment}}.emb", segment=POS_ENCODINGS),
        idx_config=f"{config.out_dir}svs_index/{{segment}}.emb_config", segment=POS_ENCODINGS),
        idx_data=f"{config.out_dir}svs_index/{{segment}}.emb_data", segment=POS_ENCODINGS),
        idx_graph=f"{config.out_dir}svs_index/{{segment}}.emb_data", segment=POS_ENCODINGS)
        #f"{config.out_dir}svs_index/idx.done"
	#f"{config.out_dir}log/faiss_build.log",

# 1.0 split all_embeddings.txt into segment embeddings for input to indexing steps
rule split_embeddings:
    input:
        f"{config.out_dir}embeddings/all.embeddings.txt"
    output:
        embeddings=expand(f"{config.out_dir}embeddings/{{segment}}.emb", segment=POS_ENCODINGS),
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
        embeddings=xpand(f"{config.out_dir}embeddings/{{segment}}.emb", segment=POS_ENCODINGS)
    output:
        idx_config=f"{config.out_dir}svs_index/{{segment}}.config", segment=POS_ENCODINGS),
        idx_data=f"{config.out_dir}svs_index/{{segment}}.data", segment=POS_ENCODINGS),
        idx_graph=f"{config.out_dir}svs_index/{{segment}}.graph", segment=POS_ENCODINGS)
        #idx_done=f"{config.out_dir}svs_index/idx.done"
    message:
        "SVS-building indexes..."
    shell:
        "conda activate svs;"
        "test ! -d {config.out_dir}svs_index/ && mkdir {config.out_dir}svs_index/;" \
        "python {config.root_dir}python/scripts/svs/build_svs_index.py" \
        " --emb_dir {config.out_dir}embeddings/" \
        " --idx_dir {config.out_dir}svs_index/" \
        " --db_samples {config.out_dir}database_hap_IDs.txt" \
        " --emb_ext emb;" \
        "touch {output.idx_done};"
