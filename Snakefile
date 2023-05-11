from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/SAS/SAS_config.yaml"
config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))



rule all:
	input:
		f"{config.log_dir}sample_IDs.log",
		f"{config.log_dir}interpolated.log",
                f"{config.log_dir}segment_boundary.log",
		f"{config.log_dir}slice.log",
		f"{config.log_dir}encode.log",
		f"{config.log_dir}model.log",
		f"{config.log_dir}embeddings.log",
		#f"{config.log_dir}faiss_build.log",
		#f"{config.log_dir}clean.log",
		#f"{config.log_dir}hap_IDs.log",
		##f"{config.out_dir}faiss_ivfpqr_idx.log",
		#f"{config.out_dir}faiss_l2_search.log"

# 0.1 create a file with all sample IDs
rule get_sample_IDs:
	input:
		vcf_file=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.root_dir}sample_IDs.txt"
	log:
		sample_IDs_log=f"{config.log_dir}sample_IDs.log"
	benchmark:
        	f"{config.benchmark_dir}sample_IDs.tsv"
	message:
		"Creating a list of all sample IDs..."
	conda:
		f"{config.conda_pmed}"	
	shell:
		"bcftools query -l {input.vcf_file} > {output.sample_IDs};"
		"cp {output.sample_IDs} {config.root_dir}database_IDs.txt;" \
		"cp {output.sample_IDs} {config.root_dir}query_IDs.txt;"

# 0.2 interpolate map
# one cm for every bp in 1kg
rule interpolate_map:
	input:
		vcf_file=f"{config.vcf_file}",
		ref_map=f"{config.ref_map}",
		interpolate_map_cpp=f"{config.cpp_src_dir}interpolate_map.cpp",
	output:
		vcf_bp=f"{config.root_dir}vcf_bp.txt",
		bin=f"{config.cpp_bin_dir}interpolate-map",
		interpolated_map=f"{config.root_dir}interpolated.map"
	log:		
		interpolated_log=f"{config.log_dir}interpolated.log"
	benchmark:
        	f"{config.benchmark_dir}interpolate.tsv"
	message:
		"Interpolating map file..."
	conda:
		f"{config.conda_pmed}"	
	shell:
		"bcftools query -f '%CHROM %POS\n' {input.vcf_file} > {output.vcf_bp};"
		"g++" \
                " {input.interpolate_map_cpp}" \
                " -I {config.cpp_include_dir}" \
                " -o {output.bin};"
		" {output.bin} {output.vcf_bp} {input.ref_map} {output.interpolated_map} > {log.interpolated_log};"

#0.3-A write segment boundary file (compile)
rule segment_boundary_file_compile:
	input:
		interpolated_map=f"{config.root_dir}interpolated.map",
		segment_boundary_map_cpp=f"{config.cpp_src_dir}segment_boundary_map.cpp",
	output:
		bin=f"{config.cpp_bin_dir}segment-boundary"
	log:
		segment_boundary_cp_log=f"{config.log_dir}config_boundary_cp.log"
	message:
		"Compiling--write segment boundary file..."
	conda:
		f"{config.conda_pmed}"	
	shell:
		"g++" \
		" {input.segment_boundary_map_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -I {config.htslib_dir}" \
		" -lhts" \
		" -o {output.bin}"

# 0.3-B write segment boundary file (execute)
rule segment_boundary_file_execute:
	input:
		interpolated_map=f"{config.root_dir}interpolated.map",
		bin=f"{config.cpp_bin_dir}segment-boundary"
	output:
		segment_boundary_file=f"{config.root_dir}segment_boundary.map"
	log:
		segment_boundary_ex_log=f"{config.log_dir}segment_boundary.log"
	benchmark:
        	f"{config.benchmark_dir}segment_bondary.tsv"
	message:
		"Executing--write segment boundary file..."
	conda:
		f"{config.conda_pmed}"	
	shell:
		"{input.bin} {input.interpolated_map} {output.segment_boundary_file} > {log.segment_boundary_ex_log}"

# 1.3 slice VCF into segments
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
		"echo 2. ---ENCODING VCF SEGMENTS---;" \
		"for vcf_f in {config.vcf_segments_dir}*.vcf.gz; do" \
		"	filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.vcf.*}};" \
		"	echo ... $seg_name;" \
		"	chrm_idx=${{seg_name#*chrm}};" \
		"	chrm_idx=${{chrm_idx%%.*}};" \
		#"	echo $chrm_idx $vcf_f {config.root_dir}sample_IDs.txt {config.encoding_file} {config.root_dir}interpolated.map;" \
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

# 4.0 run model 
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
		"python {config.model_dir}encode_samples.py" \
        	"	--encoder {config.model_checkpoint}" \
        	"	--output {config.embeddings_dir}embeddings.txt" \
        	"	--files {config.encodings_dir}*.gt" \
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

## 5.1 build FAISS indices
#rule faiss_build:
#        input:
#                slice_log=f"{config.log_dir}slice.log",
#                encode_log=f"{config.log_dir}encode.log",
#                model_log=f"{config.log_dir}model.log",
#		split_embeddings=f"{config.log_dir}embeddings.log"
#	log:
#		faiss_build_log=f"{config.log_dir}faiss_build.log"
#	message:
#		"FAISS-building indexes..."
#        conda:
#		f"{config.conda_pmed}"
#        shell:
#		"python {config.python_dir}faiss/build_faiss_index.py" \
#		"	--emb_dir {config.embeddings_dir}" \
#		"	--idx_dir {config.faiss_index-dir}"
