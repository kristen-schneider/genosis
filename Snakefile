from types import SimpleNamespace
configfile: "config_examples/config_GBR.yaml" # path to the config
config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_dir}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))


rule all:
	input:
		f"{config.vcf_file}",
		f"{config.samples_dir}GBR_sampleIDs.txt",
		f"{config.cpp_bin_dir}slice-vcf",
		f"{config.out_dir}slice.log",
		f"{config.cpp_bin_dir}encode-vcf",
		f"{config.out_dir}encode.log"

# 0. create a file with all sample IDs
# one line per sample ID
rule get_sample_IDs:
	input:
		vcf=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.samples_dir}GBR_sampleIDs.txt",
		sample_IDs_done=f"{config.samples_dir}GBR_sampleIDs.done"
	message:
		"Creating a list of all sample IDs from VCF file..."
	shell:
		"bcftools query -l {input.vcf} > {output.sample_IDs}"
		" && touch {output.sample_IDs_done}"

# 1.1 slice VCF into segments (compile)
rule slice_VCF_compile:
	input:
		slice_vcf_cpp=f"{config.cpp_src_dir}slice_vcf.cpp",
		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp"
	output:
		bin=f"{config.cpp_bin_dir}slice-vcf"
	message:
		"Compiling--slice vcf into segments..."
	shell:
		"g++" \
		" {input.slice_vcf_cpp}" \
		" {input.read_config_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -lhts" \
		" -o {output.bin}"
# 1.2 slice VCF into segments (execute)
rule slice_VCF_execute:
	input:
		bin=f"{config.cpp_bin_dir}slice-vcf",
		config_file=f"{config.cpp_configs_dir}sample.config"
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
		encode_vcf_cpp=f"{config.cpp_src_dir}encode_vcf.cpp",
		read_config_cpp=f"{config.cpp_src_dir}read_config.cpp",
		map_encodings_cpp=f"{config.cpp_src_dir}map_encodings.cpp",
		utils_cpp=f"{config.cpp_src_dir}utils.cpp"
	output:
		bin=f"{config.cpp_bin_dir}encode-vcf"
	message:
		"Compiling--encode vcf segments..."
	shell:
		"g++" \
		" {input.encode_vcf_cpp}" \
		" {input.read_config_cpp}" \
		" {input.map_encodings_cpp}" \
		" {input.utils_cpp}" \
		" -I {config.cpp_include_dir}" \
		" -lhts" \
		" -o {output.bin}"
# 2.2 encode vcf segments (execute)
rule encode_vcf_segments_execute:
	input:
		bin=f"{config.cpp_bin_dir}encode-vcf",
		sample_IDs=f"{config.samples_dir}GBR_sampleIDs.txt",
		config_file=f"{config.cpp_configs_dir}sample.config"
	output:
		encode_log=f"{config.out_dir}encode.log"
	message:
		"Executing--encode vcf segments..."
	shell:
		"for vcf_f in {config.out_dir}*.vcf; do" \
		"	filename=$(basename $vcf_f);" \
		"	seg_name=${{filename%.*}};" \
		"	./{input.bin} {input.config_file} {input.sample_IDs} $vcf_f {config.out_dir}${{seg_name}}.encoded > {output.encode_log};" \
		"done"
