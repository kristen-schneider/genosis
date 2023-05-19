from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/SAS/SAS_config.yaml"
configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"

config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

rule all:
	input:
		f"{config.log_dir}faiss_search.log"

# 1.0 search FAISS indices
rule faiss_build:
        input:
                slice_log=f"{config.log_dir}slice.log",
                encode_log=f"{config.log_dir}encode.log",
                model_log=f"{config.log_dir}model.log",
		split_embeddings=f"{config.log_dir}embeddings.log",
		faiss_index=f"{config.log_dir}faiss_build.log"
	log:
		faiss_search_log=f"{config.log_dir}faiss_search.log"
	benchmark:
		f"{config.benchmark_dir}faiss_search.tsv"
	message:
		"FAISS-searching indexes..."
	conda:
		f"{config.conda_faiss}"
	shell:
		"test ! -d {config.faiss_results_dir} && mkdir {config.faiss_results_dir};" \
		"python {config.python_dir}faiss/search_faiss_index.py" \
		"	--idx_dir {config.faiss_index_dir}" \
		"	--emb_dir {config.embeddings_dir}" \
		"	--k {config.k}" \
		"	--database_samples {config.database_IDs}" \
		"	--query_samples {config.query_IDs}" \
		"	--out_dir {config.faiss_results_dir}"
