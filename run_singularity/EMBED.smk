from types import SimpleNamespace
#
config = SimpleNamespace(**config)

shell.prefix("""
#source ~/.bashrc;
. /home/sdp/miniconda3/etc/profile.d/conda.sh;
conda activate pmed;
""")

import glob
from os.path import basename

encode_dir=f"{config.out_dir}encodings/"
pos_encodings=glob.glob(encode_dir + "*.pos")
pos_encodings=list(map(basename, pos_encodings))
pos_encodings=[".".join(p.split('.')[:-1]) for p in pos_encodings]
assert len(pos_encodings) > 0, "no positional encodings.."


rule all:
    input:
        f"{config.out_dir}zeros.out",
        f"{config.out_dir}sample_hap_IDs.txt",
        f"{config.out_dir}embeddings/model.out",
        f"{config.out_dir}embeddings/embeddings.out"

# 0.0 remove encodings with empty entries
rule remove_empty_encodings:
    output:
        f"{config.out_dir}zeros.out"
    message:
        "Removing positional encodings with empty entries"
    shell:
        "python {config.root_dir}python/scripts/check_pos_encodigs.py" \
        " --pos_dir {config.out_dir}encodings/" \
        " --pos_ext pos;" \
        "touch {config.out_dir}zeros.out;"

# 1.0 remove intermediate files
rule remove_intermediate_files:
    message:
        "Removing intermediate files after encoding."
    shell:
        "rm {config.out_dir}vcf_bp.txt;" \
        "rm -r {config.out_dir}*.vcf.*;" \

# 2.0 get hap_IDs for database samples
rule hap_IDs:
    output:
        sampleIDs=f"{config.out_dir}sample_hap_IDs.txt"
    message:
        "Getting haplotype IDs from encoding file..."
    shell:
        "for enc_f in {config.out_dir}encodings/*.gt; do" \
        " awk '{{print $1}}' $enc_f > {output.sampleIDs};" \
        " break;" \
        "done;" \
        "cp {output.sampleIDs} {config.out_dir}database_hap_IDs.txt;" \
        "cp {output.sampleIDs} {config.out_dir}query_hap_IDs.txt;"

# 4 run model 
rule model:
    input:
        f"{config.out_dir}zeros.out"
    output:
        model_out=f"{config.out_dir}embeddings/model.out"
    message:
        "Running model to create embedding vectors..."
    shell:
        "conda activate torch-cpu;"
        "echo 3. ---RUNNING MODEL---;" \
        "test ! -d {config.out_dir}embeddings/ && mkdir {config.out_dir}embeddings/;" \
        "python {config.root_dir}pytorch/encode_samples.py" \
            " --encoder {config.root_dir}last.ckpt" \
            " --output {config.out_dir}embeddings/all.embeddings.txt" \
        #" --gpu" \
        " --files {config.out_dir}encodings/*.pos" \
        " --batch-size {config.batch_size}" \
        " --num-workers {config.n_workers};" \
	"touch {output.model_out};"

# 5.0 split all_embeddings.txt into segment embeddings
rule split_embeddings:
    input:
        model_out=f"{config.out_dir}embeddings/model.out",
    output:
        embeddings_out=f"{config.out_dir}embeddings/embeddings.out"
    message:
        "Splitting full embedding file into segments..."
    shell:
        "conda activate pmed;"
        "python {config.root_dir}python/scripts/split_embeddings.py" \
        " --emb_dir {config.out_dir}embeddings/" \
        " --all_emb {config.out_dir}embeddings/all.embeddings.txt;" \
	"touch {output.embeddings_out};"
