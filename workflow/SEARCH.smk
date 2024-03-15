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

IDX_DIR=f"{config.out_dir}svs_index/"
IDX_SEGMENTS=glob.glob(IDX_DIR + "*.config/")
IDX_SEGMENTS=[i.split('/')[-2] for i in IDX_SEGMENTS]
IDX_SEGMENTS=[idx_seg.replace('.config', '') for idx_seg in IDX_SEGMENTS]
assert len(IDX_SEGMENTS) > 0, "no indexes.."
#print(IDX_SEGMENTS)

SAMPLES = []
with open(f"{config.out_dir}sample_IDs.txt") as f:
    for line in f:
        sample_name=line.strip()
        SAMPLES.append(sample_name)

assert len(SAMPLES) > 0, "no samples..."


rule all:
    input:
        f"{config.out_dir}search_time.log"
        #expand(f"{config.out_dir}svs_results/{{segment}}.knn", segment=IDX_SEGMENTS)

# 6 search SVS indices
rule svs_search:
    input:
        idx_dir=f"{config.out_dir}svs_index/"
    output:
        f"{config.out_dir}search_time.log"
    message:
        "SVS-searching indexes and aggregating results..."
    shell:
        """
        conda activate svs;
        test ! -d {config.out_dir}top_hits/ && mkdir {config.out_dir}top_hits/;
        time python \
         {config.root_dir}python/scripts/svs/search_svs_index_samples.py \
         --idx_dir {config.out_dir}svs_index/ \
         --emb_dir {config.out_dir}embeddings/ \
         --db_samples {config.out_dir}database_hap_IDs.txt \
         --query_samples {config.out_dir}query_hap_IDs.txt \
         --knn {config.k} \
         --out_dir {config.out_dir}top_hits/ > {config.out_dir}search_time.log;
       """
       # num_segments=$(tail -n 1 {config.out_dir}segment_boundary.map | awk '{{print $2}}');
       # echo ${{num_segments}} 
       # #time touch {config.out_dir}search_time.log;
       # time python \
       #  {config.root_dir}python/scripts/svs/search_svs_index_samples.py \
       #  --idx_dir {config.out_dir}svs_index/ \
       #  --emb_dir {config.out_dir}embeddings/ \
       #  --db_samples {config.out_dir}database_hap_IDs.txt \
       #  --query_samples {config.out_dir}query_hap_IDs.txt \
       #  --num_segments ${num_segments} \
       #  --knn {config.k} \
       #  --out_dir {config.out_dir}top_hits/ > {config.out_dir}search_time.log;
        #"""

## 6 search FAISS indices
#rule faiss_search:
#    input:
#        IDX_SEGMENTS=f"{config.out_dir}faiss_index/{{segment}}.idx"
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
#        " --idx {input.IDX_SEGMENTS}" \
#        " --emb_dir {config.out_dir}embeddings/" \
#        " --emb_ext emb" \
#        " --k {config.k}" \
#        " --database_samples {config.database_IDs}" \
#        " --query_samples {config.query_IDs}" \
#        " --out_dir {config.faiss_results_dir}" \
#        " --pedigree {config.pedigree}"

## 6 search SVS indices
#rule svs_search:
#    input:
#        idx_segments=f"{config.out_dir}svs_index/{{segment}}.config/"
#        #idx_done=f"{config.out_dir}svs_index/idx.done"
#    output:
#        knn_segments=f"{config.out_dir}svs_results/{{segment}}.knn"
#        #knn_segments=f"{config.out_dir}svs_results/{{segment}}.knn"
#    message:
#        "SVS-searching indexes..."
#    shell:
#        """
#        conda activate svs;
#        test ! -d {config.out_dir}svs_results/ && mkdir {config.out_dir}svs_results/;
#
#        {{ echo {input.idx_segments} &&
#        time \
#        python -W ignore\
#         {config.root_dir}python/scripts/svs/search_svs_index.py\
#         --seg_idx {input.idx_segments}\
#         --emb_dir {config.out_dir}embeddings/\
#         --emb_ext emb\
#         --db_samples {config.out_dir}database_hap_IDs.txt\
#         --q_samples {config.out_dir}query_hap_IDs.txt\
#         --knn {config.k}\
#         --out_dir {config.out_dir}svs_results/; }} 2>> {config.out_dir}search.log;
#        """
