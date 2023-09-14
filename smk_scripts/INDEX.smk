from types import SimpleNamespace
#
config = SimpleNamespace(**config)

shell.prefix("""
. /opt/conda/etc/profile.d/conda.sh
conda activate pmed;
""")

import glob
from os.path import basename

emb_dir=f"{config.out_dir}embeddings/"
emb_segments=glob.glob(emb_dir + "*.emb")
emb_segments=list(map(basename, emb_segments))
emb_segments=[".".join(e.split('.')[:-1]) for e in emb_segments]
assert len(emb_segments) > 0, "no embeddings.."

rule all:
    input:
        f"{config.out_dir}svs_index/idx.done"
        #expand(f"{config.out_dir}faiss_index/{{segment}}.idx", segment=emb_segments),
        #expand(f"{config.out_dir}svs_index/{{segment}}.emb_config/", segment=emb_segments),
	#f"{config.out_dir}log/faiss_build.log",

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
        emb_dir=f"{config.out_dir}embeddings/"
        #emb_segments=f"{config.out_dir}svs_index/{{segment}}.emb"
    output:
        idx_done=f"{config.out_dir}svs_index/idx.done"
        #idx_segments=f"{config.out_dir}svs_index/{{segment}}.emb_config/"
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
