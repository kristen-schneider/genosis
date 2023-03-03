from glob import glob
from os.path import basename, splitext
from types import SimpleNamespace

#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/simulated/config_snakemake.yaml"
config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_dir}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

rule all:
	input:
		f"{config.vcf_file}",
		f"{config.data_dir}samples_IDs.txt",
		f"{config.data_dir}samples.log",
		f"{config.data_dir}vcf.bps",
                f"{config.data_dir}interpolated.map",
                f"{config.data_dir}interpolated.log",
		f"{config.data_dir}segment_boundary.log",
                f"{config.data_dir}segment_boundary.map",
		f"{config.out_dir}slice.log",
		f"{config.out_dir}encode.log",
		f"{config.data_dir}samples_hap_IDs.txt",
		f"{config.data_dir}database_hap_IDs.txt",
		f"{config.data_dir}query_hap_IDs.txt",
		f"{config.cpp_bin_dir}faiss-l2-build",
		f"{config.out_dir}faiss_idx.log",
		f"{config.cpp_bin_dir}faiss-l2-search",
		f"{config.out_dir}faiss_search.log"		

# 0.1 create a file with all sample IDs
# one line per sample ID
rule get_sample_IDs:
	input:
		vcf=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.data_dir}samples_IDs.txt",
		sample_IDs_done=f"{config.data_dir}samples.log"
	message:
		"Creating a list of all sample IDs from VCF file..."
	conda:
		"{config.conda_dir}"	
	shell:
		"bcftools query -l {input.vcf} > {output.sample_IDs}"
		" && touch {output.sample_IDs_done};" \
		"cp {output.sample_IDs} {config.data_dir}database_IDs.txt;" \
		"cp {output.sample_IDs} {config.data_dir}query_IDs.txt;"

# 0.2 interpolate map
# one cm for every bp in 1kg
rule interpolate_map:
	input:
		vcf_file=f"{config.vcf_file}",
		ref_map=f"{config.ref_map}",
		vcf_bp_py=f"{config.python_dir}utils/vcf_bps.py",
		interpolate_py=f"{config.python_dir}utils/interpolate_map_1kg.py"
	output:
		vcf_bps=f"{config.data_dir}vcf.bps",
		interpolated_map=f"{config.data_dir}interpolated.map",
		interpolated_log=f"{config.data_dir}interpolated.log",
	message:
		"Interpolating map file..."
	conda:
		"{config.conda_dir}"	
	shell:
		"python {input.vcf_bp_py} {input.vcf_file} {output.vcf_bps};"
		"python {input.interpolate_py} {output.vcf_bps} {input.ref_map} {output.interpolated_map};"
		"touch {output.interpolated_log};"

# 1.1 write segment boundary file (compile)
rule segment_boundary_file_compile:
	input:
		interpolated_map=f"{config.data_dir}interpolated.map",
                interpolated_log=f"{config.data_dir}interpolated.log",
		vcf_file=f"{config.vcf_file}",
		main_segment_cpp=f"{config.cpp_src_dir}main_segment.cpp",
		segment_boundary_map_cpp=f"{config.cpp_src_dir}segment_boundary_map.cpp",
		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
		read_map_cpp=f"{config.cpp_src_dir}read_map.cpp",
		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}segment-boundary"
	message:
		"Compiling--write segment boundary file..."
	conda:
		"{config.conda_dir}"	
	shell:
		"g++" \
		" {input.main_segment_cpp}" \
		" {input.segment_boundary_map_cpp}" \
		" {input.read_config_cpp}" \
		" {input.read_map_cpp}" \
		" {input.utils_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -I {config.htslib_dir}" \
		" -lhts" \
		" -o {output.bin}"

# 1.2 write segment boundary file (execute)
rule segment_boundary_file_execute:
	input:
		bin=f"{config.cpp_bin_dir}segment-boundary",
		config_file=f"{config.config_file}"
	output:
		segment_boundary_log=f"{config.data_dir}segment_boundary.log",
		segment_boundary_file=f"{config.data_dir}segment_boundary.map"
	message:
		"Executing--write segment boundary file..."
	conda:
		"{config.conda_dir}"	
	shell:
		"{input.bin} {input.config_file} > {output.segment_boundary_log}"

# 1.3 slice VCF into segments
rule slice_VCF:
	input:
		vcf_file=f"{config.vcf_file}",
		segment_boundary_file=f"{config.data_dir}segment_boundary.map",
		segment_boundary_log=f"{config.data_dir}segment_boundary.log",
	output:
		slice_log=f"{config.out_dir}slice.log"
	message:
		"Slicing VCF into segments..."
	conda:
		"{config.conda_dir}"	
	shell:
		"echo 1. ---SLICING VCF INTO SEGMENTS---;" \
		"while IFs= read -r segment start_bp end_bp; do" \
		"	echo slicing segment ${{segment}} >> {output.slice_log};" \
		"	bcftools view -h {input.vcf_file} > {config.out_dir}segment.${{segment}}.vcf;" \
		"	tabix {input.vcf_file} chr8:${{start_bp}}-${{end_bp}} >> {config.out_dir}segment.${{segment}}.vcf;" \
		" done < {input.segment_boundary_file};" \

# 2.1 encode genoypes for VCF segments (compile)
rule encode_compile:
	input:
		slice_log=f"{config.out_dir}slice.log",
		main_encode_cpp=f"{config.cpp_src_dir}main_encode.cpp",
                encode_segment_cpp=f"{config.cpp_src_dir}encode_segment.cpp",
                read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
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
		" {input.read_config_cpp}" \
		" {input.read_map_cpp}" \
		" {input.map_encodings_cpp}" \
		" {input.utils_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -I {config.htslib_dir}" \
		" -lhts" \
		" -o {output.bin}"
		

# 2.2 encode genotypes for VCF segments (execute)
rule encode_execute:
	input:
		bin=f"{config.cpp_bin_dir}encode",
		config_file=f"{config.config_file}",
	output:
		encode_log=f"{config.out_dir}encode.log"
	message:
                "Compiling--encoding segments..."
	conda:
		"{config.conda_dir}"	
	shell:
		"echo 2. ---ENCODING VCF SEGMENTS---;" \
		"bash {config.bash_dir}encode.sh" \
                " ./{input.bin} {input.config_file} {config.out_dir}" 
		" > {config.out_dir}encode.log"

# 3.0 get hap_IDs for database samples
rule hap_IDs:
	input:
		encode_log=f"{config.out_dir}encode.log"
	output:
		sample_hap_ids=f"{config.data_dir}samples_hap_IDs.txt",
		database_hap_ids=f"{config.data_dir}database_hap_IDs.txt",
		query_hap_ids=f"{config.data_dir}query_hap_IDs.txt"
	message:
		"Getting haplotype IDs from encoding file..."
	conda:
		"{config.conda_dir}"	
	shell:
		"for enc_f in {config.out_dir}*.gt; do" \
		"	awk '{{print $1}}' $enc_f > {output.sample_hap_ids};" \
		"	break;" \
		"done;" \
		"cp {output.sample_hap_ids} {output.database_hap_ids};" \
		"cp {output.sample_hap_ids} {output.query_hap_ids};"

# 3.1 build faiss index for encoding segments (compile)
rule build_faiss_index_compile:
	input:
		encode_log=f"{config.out_dir}encode.log",
		faiss_l2_build_cpp=f"{config.cpp_src_dir}faiss_l2_build.cpp",
		faiss_utils_cpp=f"{config.cpp_src_dir}faiss_utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}faiss-l2-build"
	message:
		"Compiling--building faiss indices for all segments..."
	conda:
		"{config.conda_dir}"	
	shell:
		"g++" \
		" {input.faiss_l2_build_cpp}" \
		" {input.faiss_utils_cpp}" \
		" -I {config.conda_dir}/include/" \
		" -I {config.cpp_include_dir}" \
		" -L {config.conda_dir}/lib/" \
		" -lfaiss" \
		" -o {output.bin}"
# 3.2 build faiss index for encoding segments (execute)
rule build_faiss_index_execute:
        input:
                bin=f"{config.cpp_bin_dir}faiss-l2-build",
                database_hap_IDs=f"{config.data_dir}samples_hap_IDs.txt"
        output:
                faiss_idx_log=f"{config.out_dir}faiss_idx.log"
        message:
                "Executing--building faiss indices for all segments..."
	conda:
		"{config.conda_dir}"	
	shell:
		"for enc_f in {config.out_dir}*.gt; do" \
                "       filename=$(basename $enc_f);" \
                "       seg_name=${{filename%.*}};" \
                "       echo SEGMENT: $seg_name;" \
                "       ./{input.bin} {input.database_hap_IDs} $enc_f {config.out_dir}${{seg_name}}.faiss.idx" \
                "        >> {output.faiss_idx_log};" \
		"done;"

# 4.1 search faiss index for encodig segments (compile)
rule search_faiss_index_compile:
	input:
		faiss_idx_log=f"{config.out_dir}faiss_idx.log",
		faiss_l2_search_cpp=f"{config.cpp_src_dir}faiss_l2_search.cpp",
		faiss_utils_cpp=f"{config.cpp_src_dir}faiss_utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}faiss-l2-search"
	message:
		"Compiling--searching faiss indices for all segments..."
	conda:
		"{config.conda_dir}"	
	shell:
		"g++" \
		" {input.faiss_l2_search_cpp}" \
		" {input.faiss_utils_cpp}" \
		" -I {config.conda_dir}/include/" \
		" -I {config.cpp_include_dir}" \
		" -L {config.conda_dir}/lib/" \
		" -lfaiss" \
		" -o {output.bin}"
# 4.2 search faiss index for encoding segments (execute)
rule search_faiss_index_execute:
        input:
                bin=f"{config.cpp_bin_dir}faiss-l2-search",
                database_hap_IDs=f"{config.data_dir}samples_hap_IDs.txt",
		query_hap_IDs=f"{config.data_dir}query_hap_IDs.txt",
        output:
                faiss_search_log=f"{config.out_dir}faiss_search.log"
        message:
                "Executing--searching faiss indices for all segments..."
	conda:
		"{config.conda_dir}"
	shell:
                "for idx_f in {config.out_dir}*.idx; do" \
                "       filename=$(basename $idx_f);" \
                "       seg_faiss=${{filename%.*}};" \
		"	seg_name=${{seg_faiss%.*}};" \
                "       echo SEGMENT: $seg_name;" \
		"	encoding_f=$seg_name.gt;" \
		"	search_time_start=`date +%s.%N`;" \
		"	echo START: $search_time_start;" \
		"       ./{input.bin} $idx_f {input.database_hap_IDs} {input.query_hap_IDs} {config.out_dir}${{encoding_f}} {config.k} {config.out_dir}$seg_name.faiss.out " \
               	"        >> {output.faiss_search_log};" \
		"	search_time_end=`date +%s.%N`;" \
		"	echo END: $search_time_end;" \
		#"	single_search_time=$( echo '$search_time_end - $search_time_start' | bc -l );" \
		#"	echo SINGLE SEARCH: $single_search_time;" \
		"done;" \
		#"full_end=`date +%s.%N`;" \
		#"full_time=$( '$full_start - $full_end' | bc -l );" \
		#"echo FULL TIME: $full_time;"


