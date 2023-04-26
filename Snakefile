from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/SAS/SAS_config.yaml"
config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_dir}/lib"
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
		#f"{config.log_dir}clean.log",
		#f"{config.log_dir}hap_IDs.log",
		#f"{config.cpp_bin_dir}faiss-l2-build",
		#f"{config.out_dir}faiss_l2_idx.log",
		##f"{config.cpp_bin_dir}faiss-ivfpqr-build",
		##f"{config.out_dir}faiss_ivfpqr_idx.log",
		#f"{config.cpp_bin_dir}faiss-search",
		#f"{config.out_dir}faiss_l2_search.log"
		#f"{config.out_dir}faiss_ivfpqr_search.log"

# 0.1 create a file with all sample IDs
rule get_sample_IDs:
	input:
		vcf_file=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.data_dir}sample_IDs.txt"
	log:
		sample_IDs_log=f"{config.log_dir}sample_IDs.log"
	benchmark:
        	f"{config.benchmark_dir}sample_IDs.tsv"
	message:
		"Creating a list of all sample IDs..."
	conda:
		"{config.conda_dir}"	
	shell:
		"bcftools query -l {input.vcf_file} > {output.sample_IDs};"
		"cp {output.sample_IDs} {config.data_dir}database_IDs.txt;" \
		"cp {output.sample_IDs} {config.data_dir}query_IDs.txt;"

# 0.2 interpolate map
# one cm for every bp in 1kg
rule interpolate_map:
	input:
		vcf_file=f"{config.vcf_file}",
		ref_map=f"{config.ref_map}",
		interpolate_map_cpp=f"{config.cpp_src_dir}interpolate_map.cpp",
	output:
		vcf_bp=f"{config.data_dir}vcf_bp.txt",
		bin=f"{config.cpp_bin_dir}interpolate-map",
		interpolated_map=f"{config.data_dir}interpolated.map"
	log:		
		interpolated_log=f"{config.log_dir}interpolated.log"
	benchmark:
        	f"{config.benchmark_dir}interpolate.tsv"
	message:
		"Interpolating map file..."
	conda:
		"{config.conda_dir}"	
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
		interpolated_map=f"{config.data_dir}interpolated.map",
		segment_boundary_map_cpp=f"{config.cpp_src_dir}segment_boundary_map.cpp",
	output:
		bin=f"{config.cpp_bin_dir}segment-boundary"
	log:
		segment_boundary_cp_log=f"{config.log_dir}config_boundary_cp.log"
	message:
		"Compiling--write segment boundary file..."
	conda:
		"{config.conda_dir}"	
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
		interpolated_map=f"{config.data_dir}interpolated.map",
		bin=f"{config.cpp_bin_dir}segment-boundary"
	output:
		segment_boundary_file=f"{config.data_dir}segment_boundary.map"
	log:
		segment_boundary_ex_log=f"{config.log_dir}segment_boundary.log"
	benchmark:
        	f"{config.benchmark_dir}.bondary.tsv"
	message:
		"Executing--write segment boundary file..."
	conda:
		"{config.conda_dir}"	
	shell:
		"{input.bin} {input.interpolated_map} {output.segment_boundary_file} > {log.segment_boundary_ex_log}"

# 1.3 slice VCF into segments
rule slice_VCF:
	input:
		vcf_file=f"{config.vcf_file}",
		segment_boundary_file=f"{config.data_dir}segment_boundary.map",
	log:
		slice_log=f"{config.log_dir}slice.log"
	benchmark:
        	f"{config.benchmark_dir}slice.tsv"
	message:
		"Slicing VCF into segments..."
	conda:
		"{config.conda_dir}"	
	shell:
		"echo 1. ---SLICING VCF INTO SEGMENTS---;" \
		"while IFs= read -r chrm segment start_bp end_bp; do" \
		"	echo slicing segment ${{segment}} >> {log.slice_log};" \
		"	bcftools view -h {input.vcf_file} > {config.out_dir}chrm${{chrm}}.segment.${{segment}}.vcf;" \
		"	tabix {input.vcf_file} chr${{chrm}}:${{start_bp}}-${{end_bp}} >> {config.out_dir}chrm${{chrm}}.segment.${{segment}}.vcf;" \
		"	bgzip {config.out_dir}chrm${{chrm}}.segment.${{segment}}.vcf;" \
		"	tabix -p vcf {config.out_dir}chrm${{chrm}}.segment.${{segment}}.vcf.gz;" \
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
		"{config.conda_dir}"	
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
		"{config.conda_dir}"	
	shell:
		"echo 2. ---ENCODING VCF SEGMENTS---;" \
		"for vcf_f in {config.out_dir}*.vcf.gz; do" \
		"	filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.vcf.*}};" \
		"	echo ... $seg_name;" \
		"	chrm_idx=${{seg_name#*chrm}};" \
		"	chrm_idx=${{chrm_idx%%.*}};" \
		#"	echo $chrm_idx $vcf_f {config.data_dir}sample_IDs.txt {config.encoding_file} {config.data_dir}interpolated.map;" \
		"	{input.bin} " \
		"		$chrm_idx " \
		"		$vcf_f " \
		"		{config.data_dir}sample_IDs.txt " \
		"		{config.encoding_file} " \
		"		{config.data_dir}interpolated.map " \
		"		{config.out_dir}${{seg_name}}.gt " \
		"		{config.out_dir}${{seg_name}}.pos " \
		"		{config.out_dir}${{seg_name}}.af " \
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
		"{config.conda_dir}"
	shell:
		"rm {config.data_dir}vcf_bp.txt;" \
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
		"{config.conda_dir}"	
	shell:
		"for enc_f in {config.out_dir}*.gt; do" \
		"	awk '{{print $1}}' $enc_f > {config.data_dir}sample_hap_IDs.txt;" \
		"	break;" \
		"done;" \
		"cp {config.data_dir}sample_hap_IDs.txt {config.data_dir}database_hap_IDs.txt;" \
		"cp {config.data_dir}sample_hap_IDs.txt {config.data_dir}query_hap_IDs.txt;"

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
		"{config.model_conda_dir}"
	shell:
		"python {config.model_dir}encode_samples.py" \
        	"	--encoder {config.model_checkpoint}" \
        	"	--output {config.model_out_dir}" \
        	"	--files {config.out_dir}*.gt" \
        	"	--batch-size {config.batch_size}" \
        	"	--num-workers {config.n_workers}"
## 4.1 build faiss index (l2) for encoding segments (compile)
#rule build_l2_faiss_index_compile:
#	input:
#		encode_log=f"{config.out_dir}encode.log",
#		faiss_l2_build_cpp=f"{config.cpp_src_dir}faiss_l2_build.cpp",
#		faiss_utils_cpp=f"{config.cpp_src_dir}faiss_utils.cpp"
#	output:
#		bin=f"{config.cpp_bin_dir}faiss-l2-build"
#	message:
#		"Compiling--building faiss l2 indices for all segments..."
#	conda:
#		"{config.conda_dir}"	
#	shell:
#		"g++" \
#		" {input.faiss_l2_build_cpp}" \
#		" {input.faiss_utils_cpp}" \
#		" -I {config.conda_dir}/include/" \
#		" -I {config.cpp_include_dir}" \
#		" -L {config.conda_dir}/lib/" \
#		" -lfaiss" \
#		" -o {output.bin}"
## 4.2 build faiss index (l2) for encoding segments (execute)
#rule build_l2_faiss_index_execute:
#        input:
#                bin=f"{config.cpp_bin_dir}faiss-l2-build",
#                database_hap_IDs=f"{config.data_dir}samples_hap_IDs.txt"
#        output:
#                faiss_idx_log=f"{config.out_dir}faiss_l2_idx.log"
#        message:
#                "Executing--building faiss l2 indices for all segments..."
#	conda:
#		"{config.conda_dir}"	
#	shell:
#		"for enc_f in {config.out_dir}*.gt; do" \
#                "       filename=$(basename $enc_f);" \
#                "       seg_name=${{filename%.*}};" \
#                "       echo SEGMENT: $seg_name;" \
#                "       ./{input.bin} {input.database_hap_IDs} $enc_f {config.out_dir}${{seg_name}}.faissl2.idx" \
#                "        >> {output.faiss_idx_log};" \
#		"done;"
## 4.3 build faiss index (hnsw) for encoding segments (compile)
#rule build_hnsw_faiss_index_compile:
#        input:
#                encode_log=f"{config.out_dir}encode.log",
#                faiss_hnsw_build_cpp=f"{config.cpp_src_dir}faiss_hnsw_build.cpp",
#                faiss_utils_cpp=f"{config.cpp_src_dir}faiss_utils.cpp"
#        output:
#                bin=f"{config.cpp_bin_dir}faiss-hnsw-build"
#        message:
#                "Compiling--building faiss hnsw indices for all segments..."
#        conda:
#                "{config.conda_dir}"
#        shell:
#                "g++" \
#                " {input.faiss_hnsw_build_cpp}" \
#                " {input.faiss_utils_cpp}" \
#                " -I {config.conda_dir}/include/" \
#                " -I {config.cpp_include_dir}" \
#                " -L {config.conda_dir}/lib/" \
#                " -lfaiss" \
#                " -o {output.bin}"
## 4.4 build faiss index (hnsw) for encoding segments (execute)
#rule build_hnsw_faiss_index_execute:
#        input:
#                bin=f"{config.cpp_bin_dir}faiss-hnsw-build",
#                database_hap_IDs=f"{config.data_dir}samples_hap_IDs.txt"
#        output:
#                faiss_idx_log=f"{config.out_dir}faiss_hnsw_idx.log"
#        message:
#                "Executing--building faiss hnsw indices for all segments..."
#        conda:
#                "{config.conda_dir}"
#        shell:
#                "for enc_f in {config.out_dir}*.gt; do" \
#                "       filename=$(basename $enc_f);" \
#                "       seg_name=${{filename%.*}};" \
#                "       echo SEGMENT: $seg_name;" \
#                "       ./{input.bin} {input.database_hap_IDs} $enc_f {config.out_dir}${{seg_name}}.faisshnsw.idx" \
#                "        >> {output.faiss_idx_log};" \
#                "done;"
## 4.5 build faiss index (ivfpqr) for encoding segments (compile)
#rule build_ivfpqr_faiss_index_compile:
#        input:
#                encode_log=f"{config.out_dir}encode.log",
#                faiss_ivfpqr_build_cpp=f"{config.cpp_src_dir}faiss_ivfpqr_build.cpp",
#                faiss_utils_cpp=f"{config.cpp_src_dir}faiss_utils.cpp"
#        output:
#                bin=f"{config.cpp_bin_dir}faiss-ivfpqr-build"
#        message:
#                "Compiling--building faiss ivfpqr indices for all segments..."
#        conda:
#                "{config.conda_dir}"
#        shell:
#                "g++" \
#                " {input.faiss_ivfpqr_build_cpp}" \
#                " {input.faiss_utils_cpp}" \
#                " -I {config.conda_dir}/include/" \
#                " -I {config.cpp_include_dir}" \
#                " -L {config.conda_dir}/lib/" \
#                " -lfaiss" \
#                " -o {output.bin}"
## 4.6 build faiss index (ivfpqr) for encoding segments (execute)
#rule build_ivfpqr_faiss_index_execute:
#        input:
#                bin=f"{config.cpp_bin_dir}faiss-ivfpqr-build",
#                database_hap_IDs=f"{config.data_dir}samples_hap_IDs.txt"
#        output:
#                faiss_idx_log=f"{config.out_dir}faiss_ivfpqr_idx.log"
#        message:
#                "Executing--building faiss ivfpqr indices for all segments..."
#        conda:
#                "{config.conda_dir}"
#        shell:
#                "for enc_f in {config.out_dir}*.gt; do" \
#                "       filename=$(basename $enc_f);" \
#                "       seg_name=${{filename%.*}};" \
#                "       echo SEGMENT: $seg_name;" \
#                "       ./{input.bin} {input.database_hap_IDs} $enc_f {config.out_dir}${{seg_name}}.faissivfpqr.idx" \
#                "        >> {output.faiss_idx_log};" \
#                "done;"
#
#
#
#
## 5.1 search faiss index for encodig segments (compile)
#rule search_l2_faiss_index_compile:
#	input:
#		faiss_idx_log=f"{config.out_dir}faiss_l2_idx.log",
#		faiss_l2_search_cpp=f"{config.cpp_src_dir}faiss_l2_search.cpp",
#		faiss_utils_cpp=f"{config.cpp_src_dir}faiss_utils.cpp"
#	output:
#		bin=f"{config.cpp_bin_dir}faiss-search"
#	message:
#		"Compiling--searching faiss indices for all segments..."
#	conda:
#		"{config.conda_dir}"	
#	shell:
#		"g++" \
#		" {input.faiss_l2_search_cpp}" \
#		" {input.faiss_utils_cpp}" \
#		" -I {config.conda_dir}/include/" \
#		" -I {config.cpp_include_dir}" \
#		" -L {config.conda_dir}/lib/" \
#		" -lfaiss" \
#		" -o {output.bin}"
## 5.2 search faiss index for encoding segments (execute)
#rule search_l2_faiss_index_execute:
#        input:
#                bin=f"{config.cpp_bin_dir}faiss-search",
#                database_hap_IDs=f"{config.data_dir}samples_hap_IDs.txt",
#		query_hap_IDs=f"{config.data_dir}query_hap_IDs.txt",
#        output:
#                faiss_search_log=f"{config.out_dir}faiss_l2_search.log"
#        message:
#                "Executing--searching L2 faiss indices for all segments..."
#	conda:
#		"{config.conda_dir}"
#	shell:
#                "for idx_f in {config.out_dir}*.faissl2.idx; do" \
#                "       filename=$(basename $idx_f);" \
#                "       seg_faiss=${{filename%.*}};" \
#		"	seg_name=${{seg_faiss%.*}};" \
#                "       echo SEGMENT: $seg_name;" \
#		"	encoding_f=$seg_name.gt;" \
#		"	search_time_start=`date +%s.%N`;" \
#		"	echo START: $search_time_start;" \
#		"       ./{input.bin} $idx_f {input.database_hap_IDs} {input.query_hap_IDs} {config.out_dir}${{encoding_f}} {config.k} {config.out_dir}$seg_name.faissl2.out " \
#               	"        >> {output.faiss_search_log};" \
#		"	search_time_end=`date +%s.%N`;" \
#		"	echo END: $search_time_end;" \
#		#"	single_search_time=$( echo '$search_time_end - $search_time_start' | bc -l );" \
#		#"	echo SINGLE SEARCH: $single_search_time;" \
#		"done;" \
#		#"full_end=`date +%s.%N`;" \
#		#"full_time=$( '$full_start - $full_end' | bc -l );" \
#		#"echo FULL TIME: $full_time;"
#
