from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/AFR/AFR_config.yaml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"
configfile: "/Users/krsc0813/chr10/config_fiji.yaml"

config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

rule all:
    input:
	f"{config.log_dir}faiss_build.log",

## 5.1 build FAISS indices
#rule faiss_build:
#        input:
#                slice_log=f"{config.log_dir}slice.log",
#                encode_log=f"{config.log_dir}encode.log",
#                model_log=f"{config.log_dir}model.log",
#		split_embeddings=f"{config.log_dir}embeddings.log"
#	log:
#		faiss_build_log=f"{config.log_dir}faiss_build.log"
#	benchmark:
#		f"{config.benchmark_dir}faiss_index.tsv"
#	message:
#		"FAISS-building indexes..."
#	conda:
#		f"{config.conda_faiss}"
#	shell:
#		"echo 4. ---CREATING FAISS INDEX---;" \
#		"test ! -d {config.faiss_index_dir} && mkdir {config.faiss_index_dir};" \
#		"python {config.python_dir}faiss/build_faiss_index.py" \
#		"	--emb_dir {config.embeddings_dir}" \
#		"	--idx_dir {config.faiss_index_dir}" \
#		"	--db_samples {config.database_IDs}"
