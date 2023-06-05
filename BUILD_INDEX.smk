from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/AFR/AFR_config.yaml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"
configfile: "/Users/krsc0813/chr10_12/config_snakemake.yaml"

config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))



rule all:
	input:
		f"{config.log_dir}slice.log",
		f"{config.log_dir}encode.log",
		f"{config.log_dir}model.log",
		f"{config.log_dir}embeddings.log",
		f"{config.log_dir}faiss_build.log",

# 1 slice VCF into segments
rule slice_VCF:
	input:
		vcf_file=f"{config.vcf_file}",
		segment_boundary_file=f"{config.root_dir}segment_boundary.map",
	log:
		slice_log=f"{config.log_dir}slice.log"
	benchmark:
        	f"{config.benchmark_dir}slice.tsv"
	message:
		"Slicing VCF into segments..."
	conda:
		f"{config.conda_pmed}"	
	shell:
		"test ! -d {config.vcf_segments_dir} && mkdir {config.vcf_segments_dir};" \
		"echo 1. ---SLICING VCF INTO SEGMENTS---;" \
		"while IFs= read -r chrm segment start_bp end_bp; do" \
		"	echo slicing segment ${{segment}} >> {log.slice_log};" \
		"	bcftools view -h {input.vcf_file} > {config.vcf_segments_dir}chrm${{chrm}}.segment.${{segment}}.vcf;" \
		"	tabix {input.vcf_file} chr${{chrm}}:${{start_bp}}-${{end_bp}} >> {config.vcf_segments_dir}chrm${{chrm}}.segment.${{segment}}.vcf;" \
		"	bgzip {config.vcf_segments_dir}chrm${{chrm}}.segment.${{segment}}.vcf;" \
		"	tabix -p vcf {config.vcf_segments_dir}chrm${{chrm}}.segment.${{segment}}.vcf.gz;" \
		" done < {input.segment_boundary_file};" \

# 2.1 encode genotypes for VCF segments (compile)
rule encode_compile:
	input:
		slice_log=f"{config.log_dir}slice.log",
		main_encode_cpp=f"{config.cpp_src_dir}main_encode.cpp",
                encode_segment_cpp=f"{config.cpp_src_dir}encode_segment.cpp",
                read_map_cpp=f"{config.cpp_src_dir}read_map.cpp",
                map_encodings_cpp=f"{config.cpp_src_dir}map_encodings.cpp",
		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}encode"
	message:
		"Compiling--encoding segments..."
	conda:
		f"{config.conda_pmed}"	
	shell:
		"g++" \
		" {input.main_encode_cpp}" \
		" {input.encode_segment_cpp}" \
		" {input.read_map_cpp}" \
		" {input.map_encodings_cpp}" \
		" {input.utils_cpp} " \
		" -I {config.cpp_include_dir}" \
		" -I {config.htslib_dir}" \
		" -lhts" \
		" -o {output.bin}"
		
# 2.2 encode genotypes for VCF segments (execute)
rule encode_execute:
	input:
		bin=f"{config.cpp_bin_dir}encode",
	output:
		encode_log=f"{config.log_dir}encode.log"
	benchmark:
        	f"{config.benchmark_dir}encode.tsv"
	message:
                "Executing--encoding segments..."
	conda:
		f"{config.conda_pmed}"	
	shell:
		"test ! -d {config.encodings_dir} && mkdir {config.encodings_dir};" \
		"echo 2. ---ENCODING VCF SEGMENTS---;" \
		"for vcf_f in {config.vcf_segments_dir}*.vcf.gz; do" \
		"	filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.vcf.*}};" \
		"	chrm_idx=${{seg_name#*chrm}};" \
		"	chrm_idx=${{chrm_idx%%.*}};" \
		"	{input.bin} " \
		"		$chrm_idx " \
		"		$vcf_f " \
		"		{config.root_dir}sample_IDs.txt " \
		"		{config.encoding_file} " \
		"		{config.root_dir}interpolated.map " \
		"		{config.encodings_dir}${{seg_name}}.gt " \
		"		{config.encodings_dir}${{seg_name}}.pos " \
		"		{config.encodings_dir}${{seg_name}}.af " \
		"		 >> {output.encode_log};" \
		"done;" \

# 3.0 remove intermediate files
rule remove_intermediate_files:
	input:
		slice_log=f"{config.log_dir}slice.log",
		encode_log=f"{config.log_dir}encode.log"
	log:
		clean_log=f"{config.log_dir}clean.log"
	message:
		"Removing intermediate files after encoding."
	conda:
		f"{config.conda_pmed}"
	shell:
		"rm {config.root_dir}vcf_bp.txt;" \
		"rm -r {config.out_dir}*.vcf.*;" \


# 3.1 get hap_IDs for database samples
rule hap_IDs:
	input:
		encode_log=f"{config.log_dir}encode.log",
		clean_log=f"{config.log_dir}clean.log"
	log:
		hap_IDs_log=f"{config.log_dir}hap_IDs.log"
	message:
		"Getting haplotype IDs from encoding file..."
	conda:
		f"{config.conda_pmed}"	
	shell:
		"for enc_f in {config.encodings_dir}*.gt; do" \
		"	awk '{{print $1}}' $enc_f > {config.root_dir}sample_hap_IDs.txt;" \
		"	break;" \
		"done;" \
		"cp {config.root_dir}sample_hap_IDs.txt {config.root_dir}database_hap_IDs.txt;" \
		"cp {config.root_dir}sample_hap_IDs.txt {config.root_dir}query_hap_IDs.txt;"

# 4 run model 
rule model:
	input:
		encode_log=f"{config.log_dir}encode.log"
	log:
		model_log=f"{config.log_dir}model.log"
	benchmark:
        	f"{config.benchmark_dir}model.tsv"
	message:
		"Running model to create embedding vectors..."
	conda:
		f"{config.conda_model}"
	shell:
		"echo 3. ---RUNNING MODEL---;" \
		"test ! -d {config.embeddings_dir} && mkdir {config.embeddings_dir};" \
		"python {config.model_dir}encode_samples.py" \
        	"	--encoder {config.model_checkpoint}" \
        	"	--output {config.embeddings_dir}embeddings.txt" \
		"	--gpu" \
        	"	--files {config.encodings_dir}*.pos" \
        	"	--batch-size {config.batch_size}" \
        	"	--num-workers {config.n_workers}"

# 5.0 split all_embeddings.txt into segment embeddings
rule split_embeddings:
	input:
		slice_log=f"{config.log_dir}slice.log",
                encode_log=f"{config.log_dir}encode.log",
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
		"python {config.python_dir}faiss/split_embeddings.py" \
		"       --emb_dir {config.embeddings_dir}" \
		"	--all_emb {config.embeddings_dir}embeddings.txt"

# 5.1 build FAISS indices
rule faiss_build:
        input:
                slice_log=f"{config.log_dir}slice.log",
                encode_log=f"{config.log_dir}encode.log",
                model_log=f"{config.log_dir}model.log",
		split_embeddings=f"{config.log_dir}embeddings.log"
	log:
		faiss_build_log=f"{config.log_dir}faiss_build.log"
	benchmark:
		f"{config.benchmark_dir}faiss_index.tsv"
	message:
		"FAISS-building indexes..."
	conda:
		f"{config.conda_faiss}"
	shell:
		"echo 4. ---CREATING FAISS INDEX---;" \
		"test ! -d {config.faiss_index_dir} && mkdir {config.faiss_index_dir};" \
		"python {config.python_dir}faiss/build_faiss_index.py" \
		"	--emb_dir {config.embeddings_dir}" \
		"	--idx_dir {config.faiss_index_dir}" \
		"	--db_samples {config.database_IDs}"
