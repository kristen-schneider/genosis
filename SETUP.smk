from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/SAS/SAS_config.yaml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"

config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))



rule all:
	input:
		f"{config.log_dir}setup.log",
		f"{config.log_dir}sample_IDs.log",
		f"{config.log_dir}interpolated.log",
                f"{config.log_dir}segment_boundary.log",

# 0 setup log and benchmark dir
rule setup_log_benchmark:
	log:
		f"{config.log_dir}setup.log"	
	conda:
		f"{config.conda_pmed}"
	message:
		"Setting up log and benchmark directories..."
	shell:
		"test ! -d {config.log_dir} && mkdir {config.log_dir}" \
		"test ! -d {config.benchmark_dir} && mkdir {config.benchmark_dir}" \

# 1 create a file with all sample IDs
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

# 2 interpolate map
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

# 3-A write segment boundary file (compile)
rule segment_boundary_file_compile:
	input:
		interpolated_map=f"{config.root_dir}interpolated.map",
		segment_boundary_map_cpp=f"{config.cpp_src_dir}segment_boundary_map.cpp",
	output:
		bin=f"{config.cpp_bin_dir}segment-boundary"
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

# 3-B write segment boundary file (execute)
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
