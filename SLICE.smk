from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/SAS/SAS_config.yaml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"
configfile: "/Users/krsc0813/chr10/config_fiji.yaml"

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
        f"{config.log_dir}vcf_bp.log",
        f"{config.log_dir}segment_boundary.log",
        f"{config.log_dir}slice.log"                

# 0 setup log and benchmark dir
rule setup_log_benchmark:
    log:
        f"{config.log_dir}setup.log"	
    conda:
        f"{config.conda_pmed}"
    message:
        "Setting up log and benchmark directories..."
    shell:
        "test ! -d {config.log_dir} && mkdir {config.log_dir};" \
        "test ! -d {config.benchmark_dir} && mkdir {config.benchmark_dir};"

# 1 create a file with all sample IDs
rule get_sample_IDs:
    input:
        setup_log=f"{config.log_dir}setup.log",
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

# 2 vcf basepairs
rule vcf_bp:
    input:
        setup_log=f"{config.log_dir}setup.log",
        vcf_file=f"{config.vcf_file}"
    output:
        vcf_bp=f"{config.root_dir}vcf_bp.txt"
    log:
        vcf_bp_log=f"{config.log_dir}vcf_bp.log"
    benchmark:
        f"{config.benchmark_dir}vcf_bp.tsv"
    conda:
        f"{config.conda_pmed}"
    shell:
        "bcftools query -f '%CHROM %POS\n' {input.vcf_file} > {output.vcf_bp};"

# 3 interpolate map
# one cm for every bp in 1kg
rule interpolate_map:
    input:
        setup_log=f"{config.log_dir}setup.log",
        vcf_bp=f"{config.root_dir}vcf_bp.txt",
        ref_map=f"{config.ref_map}",
        interpolate_map_cpp=f"{config.cpp_src_dir}interpolate_map.cpp",
    output:
        bin=f"{config.cpp_bin_dir}interpolate-map"
	#interpolated_map=f"{config.root_dir}interpolated.map"
    log:
        interpolated_log=f"{config.log_dir}interpolated.log"
    benchmark:
        f"{config.benchmark_dir}interpolate.tsv"
    message:
        "Interpolating map file..."
    conda:
        f"{config.conda_pmed}"
    shell:
        "g++" \
        " {input.interpolate_map_cpp}" \
        " -I {config.cpp_include_dir}" \
        " -o {output.bin};"
        " {output.bin} {input.vcf_bp} {input.ref_map} {config.root_dir}interpolated.map > {log.interpolated_log};"

# 4-A write segment boundary file (compile)
rule segment_boundary_file_compile:
    input:
        #interpolated_map=f"{config.root_dir}interpolated.map",
        setup_log=f"{config.log_dir}setup.log",
        interpolated_log=f"{config.log_dir}interpolated.log",
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

# 4-B write segment boundary file (execute)
rule segment_boundary_file_execute:
    input:
        f"{config.log_dir}setup.log",
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
        "{input.bin} {config.root_dir}interpolated.map {output.segment_boundary_file} > {log.segment_boundary_ex_log}"

# 5 slice VCF into segments
rule slice_VCF:
    input:
        vcf_file=f"{config.vcf_file}",
        segment_boundary_file=f"{config.root_dir}segment_boundary.map"
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
        " echo segment ${{segment}} >> {log.slice_log};" \
        " bcftools view -h {input.vcf_file} > {config.vcf_segments_dir}chrm${{chrm}}.segment${{segment}}.vcf;" \
        " tabix {input.vcf_file} chr${{chrm}}:${{start_bp}}-${{end_bp}} >> {config.vcf_segments_dir}chrm${{chrm}}.segment${{segment}}.vcf;" \
        " bgzip {config.vcf_segments_dir}chrm${{chrm}}.segment${{segment}}.vcf;" \
        " tabix -p vcf {config.vcf_segments_dir}chrm${{chrm}}.segment${{segment}}.vcf.gz;" \
        " done < {input.segment_boundary_file};"
        #"while IFs= read -r chrm segment start_bp end_bp; do" \
        #"   echo slicing segment ${{segment}} >> {log.slice_log};" \
        #" done < {input.segment_boundary_file};" \
