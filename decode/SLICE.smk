from types import SimpleNamespace

#configfile: "/home/kristen/precision-medicine/example/config_singularity.yml"
#configfile: "/nfs/fs3/ToResearch/kristens/precision-medicine/example/config_decode.yml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_singularity.yml"

config = SimpleNamespace(**config)

shell.prefix("""
source ~/.bashrc;
conda activate pmed;
""")

rule all:
    input:
        f"{config.root_dir}sample_IDs.txt",
        f"{config.root_dir}vcf_bp.txt",
        f"{config.cpp_bin_dir}interpolate-map",
        f"{config.root_dir}interpolated.map",
        f"{config.root_dir}segment_boundary.map",
        f"{config.vcf_segments_dir}vcf_slices.done"

# 1 create a file with all sample IDs
rule get_sample_IDs:
	input:
		vcf_file=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.root_dir}sample_IDs.txt"
	message:
		"Creating a list of all sample IDs..."
	shell:
		"bcftools query -l {input.vcf_file} > {output.sample_IDs};"
		"cp {output.sample_IDs} {config.root_dir}database_IDs.txt;"
		"cp {output.sample_IDs} {config.root_dir}query_IDs.txt;"

# 2 vcf basepairs
rule vcf_bp:
    input:
        vcf_file=f"{config.vcf_file}"
    output:
        vcf_bp=f"{config.root_dir}vcf_bp.txt"
    shell:
        "bcftools query -f '%CHROM %POS\n' {input.vcf_file} > {output.vcf_bp};"

# 3-A interpolate map (compile)
rule interpolate_map_compile:
    input:
        interpolate_map_cpp=f"{config.cpp_src_dir}interpolate_map.cpp",
    output:
        bin=f"{config.cpp_bin_dir}interpolate-map"
    message:
        "Compiling--interpolate map file..."
    shell:
        "g++" \
        " {input.interpolate_map_cpp}" \
        " -I {config.cpp_include_dir}" \
        " -o {output.bin};"
# 3-B interpolate map (execute)
rule interpolate_map_execute:
    input:
        vcf_bp=f"{config.root_dir}vcf_bp.txt",
        ref_map=f"{config.ref_map}",
        bin=f"{config.cpp_bin_dir}interpolate-map"
    output:
        interpolated_map=f"{config.root_dir}interpolated.map"
    message:
        "Executing--interpolate map file..."
    shell:
        "{input.bin} {input.vcf_bp} {input.ref_map} {config.root_dir}interpolated.map"

# 4-A write segment boundary file (compile)
rule segment_boundary_file_compile:
    input:
        segment_boundary_map_cpp=f"{config.cpp_src_dir}segment_boundary_map.cpp",
    output:
        bin=f"{config.cpp_bin_dir}segment-boundary"
    message:
        "Compiling--write segment boundary file..."
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
        interpolated_map=f"{config.root_dir}interpolated.map",
        bin=f"{config.cpp_bin_dir}segment-boundary"
    output:
        segment_boundary_file=f"{config.root_dir}segment_boundary.map"
    message:
        "Executing--write segment boundary file..."
    shell:
        "{input.bin} {input.interpolated_map} {output.segment_boundary_file}"

# 5 slice VCF into segments
rule slice_VCF:
    input:
        vcf_file=f"{config.vcf_file}",
        segment_boundary_file=f"{config.root_dir}segment_boundary.map"
    output:
        vcf_slices_done=f"{config.vcf_segments_dir}vcf_slices.done"
    message:
        "Slicing VCF into segments..."
    shell:
        "test ! -d {config.vcf_segments_dir} && mkdir {config.vcf_segments_dir};" \
        "echo 1. ---SLICING VCF INTO SEGMENTS---;" \
        "while IFs= read -r chrm segment start_bp end_bp; do" \
        " echo segment ${{segment}};" \
        " bcftools view -h {input.vcf_file} > {config.vcf_segments_dir}chrm${{chrm}}.segment${{segment}}.vcf;" \
        " tabix {input.vcf_file} chr${{chrm}}:${{start_bp}}-${{end_bp}} >> {config.vcf_segments_dir}chrm${{chrm}}.segment${{segment}}.vcf;" \
        " bgzip {config.vcf_segments_dir}chrm${{chrm}}.segment${{segment}}.vcf;" \
        " tabix -p vcf {config.vcf_segments_dir}chrm${{chrm}}.segment${{segment}}.vcf.gz;" \
        " done < {input.segment_boundary_file};"
        " touch {output.vcf_slices_done}"
