from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
configfile: "/home/sdp/precision-medicine/example/config_singularity.yml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/AFR/AFR_config.yaml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/Users/krsc0813/chr10/config_fiji.yaml"
#configfile: "/Users/krsc0813/chr15_20/config_snakemake.yaml"
#configfile: "/Users/krsc0813/AFR_pedigree/config_AFR.yaml"

config = SimpleNamespace(**config)

#LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
#shell.prefix("""
#set -euo pipefail;
#export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
#""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

import glob
from os.path import basename

#idx_dir=f"{config.svs_index_dir}"
#idx_segments=glob.glob(idx_dir + "*.emb_config/")
#idx_segments=list(map(basename, idx_segments))
#idx_segments=[".".join(i.split('.')[:-1]) for i in idx_segments]
#assert len(idx_segments) > 0, "no indexes.."

emb_dir=f"{config.embeddings_dir}"
emb_segments=glob.glob(emb_dir + "*.emb")
emb_segments=list(map(basename, emb_segments))
emb_segments=[".".join(e.split('.')[:-1]) for e in emb_segments]
assert len(emb_segments) > 0, "no embeddings.."

rule all:
    input:
        #expand(f"{config.faiss_results_dir}{{segment}}.knn", segment=idx_segments)
        expand(f"{config.svs_results_dir}{{segment}}.knn", segment=emb_segments)

## 6 search FAISS indices
#rule faiss_search:
#    input:
#        idx_segments=f"{config.faiss_index_dir}{{segment}}.idx"
#    output:
#        knn_segments=f"{config.faiss_results_dir}{{segment}}.knn"
#    message:
#        "FAISS-searching indexes..."
#    conda:
#        f"{config.conda_faiss}"
#    shell:
#        "echo 5. ---SEARCHING FAISS INDEX---;" \
#        "test ! -d {config.faiss_results_dir} && mkdir {config.faiss_results_dir};" \
#        "python {config.python_dir}faiss/search_faiss_index.py" \
#        " --idx {input.idx_segments}" \
#        " --emb_dir {config.embeddings_dir}" \
#        " --emb_ext emb" \
#        " --k {config.k}" \
#        " --database_samples {config.database_IDs}" \
#        " --query_samples {config.query_IDs}" \
#        " --out_dir {config.faiss_results_dir}" \
#        " --pedigree {config.pedigree}"

# 6 search SVS indices
rule svs_search:
    input:
        #idx_segments=f"{config.svs_index_dir}"
        idx_done=f"{config.svs_index_dir}idx.done"
    output:
        knn_segments=f"{config.svs_results_dir}{{segment}}.knn"
    message:
        "SVS-searching indexes..."
    conda:
        f"{config.conda_base}"
    shell:
        "echo 5. --SEARCHING SVS INDEX---;" \
        "test ! -d {config.svs_results_dir} && mkdir {config.svs_results_dir};" \
        "python {config.python_dir}svs/search_svs_index.py" \
        " --idx_dir {config.svs_index_dir}" \
        " --emb_dir {config.embeddings_dir}" \
        " --emb_ext emb" \
        " --db_samples {config.database_hap_IDs}" \
        " --q_samples {config.query_hap_IDs}" \
        " --k 20" \
        " --out_dir {config.svs_results_dir}"
