from types import SimpleNamespace

#configfile: "/nfs/fs3/ToResearch/kristens/precision-medicine/example/config_singularity.yml"
configfile: "/nfs/fs3/ToResearch/kristens/precision-medicine/example/config_decode.yml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"

config = SimpleNamespace(**config)

shell.prefix("""
source ~/.bashrc;
conda activate pmed;
""")

import glob
from os.path import basename

emb_dir=f"{config.embeddings_dir}"
emb_segments=glob.glob(emb_dir + "*.emb")
emb_segments=list(map(basename, emb_segments))
emb_segments=[".".join(e.split('.')[:-1]) for e in emb_segments]
assert len(emb_segments) > 0, "no embeddings.."

rule all:
    input:
        f"{config.svs_index_dir}idx.done"
        #expand(f"{config.faiss_index_dir}{{segment}}.idx", segment=emb_segments),
        #expand(f"{config.svs_index_dir}{{segment}}.emb_config/", segment=emb_segments),
	#f"{config.log_dir}faiss_build.log",

## 5.1 build FAISS indices
#rule faiss_build:
#    input:
#        emb_segments=f"{config.embeddings_dir}{{segment}}.emb"
#    output:
#        idx_segments=f"{config.faiss_index_dir}{{segment}}.idx"
#    message:
#        "FAISS-building indexes..."
#    conda:
#        f"{config.conda_faiss}"
#    shell:
#        "echo 4. ---CREATING FAISS INDEX---;" \
#        "test ! -d {config.faiss_index_dir} && mkdir {config.faiss_index_dir};" \
#        "python {config.python_dir}faiss/build_faiss_index.py" \
#        " --emb {input.emb_segments}" \
#        " --idx_dir {config.faiss_index_dir}" \
#        " --db_samples {config.database_IDs}" \
#        " --pedigree {config.pedigree}"
#

# 5.2 build SVS indices
rule svs_build:
	input:
		emb_dir=f"{config.embeddings_dir}"
	output:
		idx_done=f"{config.svs_index_dir}idx.done"
	message:
		"SVS-building indexes..."
	shell:
		"conda activate base;"
		"echo 4. ---CREATING SVS INDEX--;" \
		"test ! -d {config.svs_index_dir} && mkdir {config.svs_index_dir};" \
		"python {config.python_dir}svs/build_svs_index.py" \
		" --emb_dir {config.embeddings_dir}" \
		" --idx_dir {config.svs_index_dir}" \
		" --db_samples {config.database_hap_IDs}" \
		" --emb_ext emb;" \
		"touch {output.idx_done};"
