from types import SimpleNamespace
configfile: "notes/ancestry_configs/config_EAS_snakemake.yaml" # path to the config
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
		f"{config.cpp_bin_dir}slice-vcf",
		f"{config.out_dir}slice.log",
		f"{config.cpp_bin_dir}encode-vcf",
		f"{config.out_dir}encode.log",
		f"{config.data_dir}plink.log",
		f"{config.out_dir}distance.log"	,
		f"{config.data_dir}aggregate.log"	
	
# 0. create a file with all sample IDs
# one line per sample ID
rule get_sample_IDs:
	input:
		vcf=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.data_dir}samples_IDs.txt",
		sample_IDs_done=f"{config.data_dir}samples.log"
	message:
		"Creating a list of all sample IDs from VCF file..."
	shell:
		"bcftools query -l {input.vcf} > {output.sample_IDs}"
		" && touch {output.sample_IDs_done}"

# 1.1 slice VCF into segments (compile)
rule slice_VCF_compile:
	input:
		main_slice_cpp=f"{config.cpp_src_dir}main_slice.cpp",
		slice_vcf_cpp=f"{config.cpp_src_dir}slice_vcf.cpp",
		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
		read_map_cpp=f"{config.cpp_src_dir}read_map.cpp",
		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}slice-vcf"
	message:
		"Compiling--slice vcf into segments..."
	shell:
		"g++ -std=c++11" \
		" {input.main_slice_cpp}" \
		" {input.slice_vcf_cpp}" \
		" {input.read_config_cpp}" \
		" {input.read_map_cpp}" \
		" {input.utils_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -I {config.htslib_dir}" \
		" -lhts" \
		" -o {output.bin}"
# 1.2 slice VCF into segments (execute)
rule slice_VCF_execute:
	input:
		bin=f"{config.cpp_bin_dir}slice-vcf",
		config_file=f"{config.config_file}"
	output:
		slice_log=f"{config.out_dir}slice.log"
	message:
		"Executing--slice vcf into segments..."
	shell:
		"./{input.bin} {input.config_file} > {output.slice_log}"

# 2.1 encode vcf segments (compile)
rule encode_vcf_segments_compile:
	input:
		slice_log=f"{config.out_dir}slice.log",
		main_encode_cpp=f"{config.cpp_src_dir}main_encode.cpp",
		encode_vcf_cpp=f"{config.cpp_src_dir}encode_vcf.cpp",
		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
		map_encodings_cpp=f"{config.cpp_src_dir}map_encodings.cpp",
		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}encode-vcf"
	message:
		"Compiling--encode vcf segments..."
	shell:
		"g++ -std=c++11" \
		" {input.main_encode_cpp}" \
		" {input.encode_vcf_cpp}" \
		" {input.read_config_cpp}" \
		" {input.map_encodings_cpp}" \
		" {input.utils_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -I {config.htslib_dir}" \
		" -lhts" \
		" -o {output.bin}"
# 2.2 encode vcf segments (execute)
rule encode_vcf_segments_execute:
	input:
		bin=f"{config.cpp_bin_dir}encode-vcf",
		sample_IDs=f"{config.data_dir}samples_IDs.txt",
		config_file=f"{config.config_file}"
	output:
		encode_log=f"{config.out_dir}encode.log"
	message:
		"Executing--encode vcf segments..."
	shell:
		"for vcf_f in {config.out_dir}*.vcf; do" \
		"	filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.*}};" \
		"	./{input.bin} {input.config_file} $vcf_f {config.out_dir}${{seg_name}}.encoded {config.out_dir}${{seg_name}}.position;" \
		"done;"
		"touch {output.encode_log};"

# 3 run plink on full vcf
rule plink:
	input:
		vcf=f"{config.vcf_file}"
	output:
		plink_done=f"{config.data_dir}plink.log"
	message:
		"Running plink --genome on full vcf file"
	shell:
		"plink --vcf {input.vcf} --genome --out {config.data_dir}plink"

# 4.1 compute euclidean distance for all segments
rule compute_segment_distance:
	input:
		encode_log=f"{config.out_dir}encode.log",
		query_file=f"{config.query_file}"
	output:
		distance_log=f"{config.out_dir}distance.log"
	message:
		"Computing Euclidean distance for query against all segments"
	shell:
		"for encoded_f in {config.out_dir}*.encoded; do" \
		"	filename=$(basename $encoded_f);" \
                "	seg_name=${{filename%.*}};" \
		"	echo $encoded_f >> {output.distance_log};" \
		"	python {config.python_dir}distance/compute_segment_distance.py --encoded_file $encoded_f --query_file {input.query_file} > {config.out_dir}${{seg_name}}.dist;" \
		"done"
# 4.2 aggregate euclidean distance for all segments
rule aggregate_segment_distance:
	input:
		distance_log=f"{config.out_dir}distance.log"
	output:
		aggregate_log=f"{config.data_dir}aggregate.log"
	message:
		"Aggregating all distance files for all queries"
	shell:
		"num_segments=$(ls {config.out_dir}*dist | wc -l)"
		" && python {config.python_dir}distance/aggregate_segment_distance.py --out_dir {config.out_dir} --ext dist --num_seg $num_segments > {config.data_dir}aggregate.txt"
		" && touch {config.data_dir}aggregate.log"
