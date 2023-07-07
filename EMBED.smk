from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/AFR/AFR_config.yaml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/Users/krsc0813/chr10/config_fiji.yaml"
#configfile: "/Users/krsc0813/chr10_12/config_snakemake.yaml"
configfile: "/Users/krsc0813/AFR_pedigree/AFR_config.yaml"

config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

import glob
from os.path import basename

encode_dir=f"{config.encodings_dir}"
pos_encodings=glob.glob(encode_dir + "*.pos")
pos_encodings=list(map(basename, pos_encodings))
pos_encodings=[".".join(p.split('.')[:-1]) for p in pos_encodings]
assert len(pos_encodings) > 0, "no positional encodings.."


rule all:
    input:
        f"{config.log_dir}zeros.log",
        f"{config.root_dir}sample_hap_IDs.txt",
        f"{config.log_dir}model.log",
        f"{config.log_dir}embeddings.log"

# 0.0 remove encodings with empty entries
rule remove_empty_encodings:
    output:
        f"{config.log_dir}zeros.log"
    message:
        "Removing positional encodings with empty entries"
    conda:
        f"{config.conda_pmed}"
    shell:
        "python {config.python_dir}check_pos_encodigs.py" \
        " --pos_dir {config.encodings_dir}" \
        " --pos_ext pos;" \
        "touch {config.log_dir}zeros.log;"

# 1.0 remove intermediate files
rule remove_intermediate_files:
    message:
        "Removing intermediate files after encoding."
    conda:
        f"{config.conda_pmed}"
    shell:
        "rm {config.root_dir}vcf_bp.txt;" \
        "rm -r {config.out_dir}*.vcf.*;" \

# 2.0 get hap_IDs for database samples
rule hap_IDs:
    output:
        sampleIDs=f"{config.root_dir}sample_hap_IDs.txt"
    message:
        "Getting haplotype IDs from encoding file..."
    conda:
        f"{config.conda_pmed}"
    shell:
        "for enc_f in {config.encodings_dir}*.gt; do" \
        " awk '{{print $1}}' $enc_f > {output.sampleIDs};" \
        " break;" \
        "done;" \
        "cp {output.sampleIDs} {config.root_dir}database_hap_IDs.txt;" \
        "cp {output.sampleIDs} {config.root_dir}query_hap_IDs.txt;"

# 4 run model 
rule model:
    input:
        f"{config.log_dir}zeros.log"
    log:
        model_log=f"{config.log_dir}model.log"
    message:
        "Running model to create embedding vectors..."
    conda:
        f"{config.conda_model}"
    shell:
        "echo 3. ---RUNNING MODEL---;" \
        "test ! -d {config.embeddings_dir} && mkdir {config.embeddings_dir};" \
	"python {config.model_dir}encode_samples.py" \
    	"	--encoder {config.model_checkpoint}" \
    	"	--output {config.embeddings_dir}chrm10.embeddings.txt" \
	"	--gpu" \
        "	--files {config.encodings_dir}*.pos" \
        "	--batch-size {config.batch_size}" \
        "	--num-workers {config.n_workers}"

# 5.0 split all_embeddings.txt into segment embeddings
rule split_embeddings:
	input:
                model_log=f"{config.log_dir}model.log",
	log:
		embeddings_log=f"{config.log_dir}embeddings.log"
	benchmark:
                f"{config.benchmark_dir}split_embeddings.tsv"
	message:
		"Splitting full embedding file into segments..."
	conda:
		f"{config.conda_pmed}"
	shell:
		"python {config.python_dir}split_embeddings.py" \
		"       --emb_dir {config.embeddings_dir}" \
		"	--all_emb {config.embeddings_dir}chrm10.embeddings.txt"

