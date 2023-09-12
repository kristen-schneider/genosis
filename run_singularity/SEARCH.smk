from types import SimpleNamespace
#
config = SimpleNamespace(**config)

shell.prefix("""
source ~/.bashrc;
conda activate pmed;
""")

import glob
from os.path import basename

#idx_dir=f"{config.out_dir}svs_index/"
#idx_segments=glob.glob(idx_dir + "*.emb_config/")
#idx_segments=list(map(basename, idx_segments))
#idx_segments=[".".join(i.split('.')[:-1]) for i in idx_segments]
#assert len(idx_segments) > 0, "no indexes.."

emb_dir=f"{config.out_dir}embeddings/"
emb_segments=glob.glob(emb_dir + "*.emb")
emb_segments=list(map(basename, emb_segments))
emb_segments=[".".join(e.split('.')[:-1]) for e in emb_segments]
assert len(emb_segments) > 0, "no embeddings.."

rule all:
    input:
        #expand(f"{config.faiss_results_dir}{{segment}}.knn", segment=idx_segments)
        expand(f"{config.out_dir}svs_results/{{segment}}.knn", segment=emb_segments)

## 6 search FAISS indices
#rule faiss_search:
#    input:
#        idx_segments=f"{config.out_dir}faiss_index/{{segment}}.idx"
#    output:
#        knn_segments=f"{config.faiss_results_dir}{{segment}}.knn"
#    message:
#        "FAISS-searching indexes..."
#    conda:
#        f"{config.conda_faiss}"
#    shell:
#        "echo 5. ---SEARCHING FAISS INDEX---;" \
#        "test ! -d {config.faiss_results_dir} && mkdir {config.faiss_results_dir};" \
#        "python {config.root_dir}python/scripts/faiss/search_faiss_index.py" \
#        " --idx {input.idx_segments}" \
#        " --emb_dir {config.out_dir}embeddings/" \
#        " --emb_ext emb" \
#        " --k {config.k}" \
#        " --database_samples {config.database_IDs}" \
#        " --query_samples {config.query_IDs}" \
#        " --out_dir {config.faiss_results_dir}" \
#        " --pedigree {config.pedigree}"

# 6 search SVS indices
rule svs_search:
    input:
        #idx_segments=f"{config.out_dir}svs_index/"
        idx_done=f"{config.out_dir}svs_index/idx.done"
    output:
        knn_segments=f"{config.out_dir}svs_results/{{segment}}.knn"
    message:
        "SVS-searching indexes..."
    shell:
	"conda activate base;"
        "echo 5. --SEARCHING SVS INDEX---;" \
        "test ! -d {config.out_dir}svs_results/ && mkdir {config.out_dir}svs_results/;" \
        "python {config.root_dir}python/scripts/svs/search_svs_index.py" \
        " --idx_dir {config.out_dir}svs_index/" \
        " --emb_dir {config.out_dir}embeddings/" \
        " --emb_ext emb" \
        " --db_samples {config.out_dir}database_hap_IDs.txt" \
        " --q_samples {config.out_dir}query_hap_IDs.txt" \
        " --k {config.k}" \
        " --out_dir {config.out_dir}svs_results/"
