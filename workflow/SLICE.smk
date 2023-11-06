from types import SimpleNamespace
config = SimpleNamespace(**config)

#shell.prefix("""
#. /opt/conda/etc/profile.d/conda.sh
#conda activate pmed;
#""")

shell.prefix("""
set -e pipefail
source  ~/.bashrc
conda activate pmed;
""")

rule all:
    input:
        f"{config.out_dir}sample_IDs.txt",
        f"{config.out_dir}vcf_segments/vcf_slices.done"

# 1.0 create a file with all sample IDs
rule get_sample_IDs:
    input:
        vcf_file=f"{config.vcf_file}"
    output:
        sample_IDs=f"{config.out_dir}sample_IDs.txt"
    message:
        "Creating a list of all sample IDs..."
    shell:
        """
        bcftools query -l {input.vcf_file} > {output.sample_IDs};
	cp {output.sample_IDs} {config.out_dir}database_IDs.txt;
	cp {output.sample_IDs} {config.out_dir}query_IDs.txt;
        """

# 2.0 vcf basepairs
rule vcf_bp:
    input:
        vcf_file=f"{config.vcf_file}"
    output:
        vcf_bp=f"{config.out_dir}vcf_bp.txt"
    shell:
        """
        bcftools query -f '%CHROM %POS\n' {input.vcf_file} > {output.vcf_bp};
        """

# 3.1 interpolate map (compile)
rule interpolate_map_compile:
    input:
        interpolate_map_cpp=f"{config.root_dir}cpp/src/interpolate_map.cpp",
    output:
        bin=f"{config.root_dir}cpp/bin/interpolate-map"
    message:
        "Compiling--interpolate map file..."
    shell:
        """
        g++\
         {input.interpolate_map_cpp}\
         -I {config.root_dir}cpp/include/\
         -o {output.bin};
        """

# 3.2 interpolate map (execute)
rule interpolate_map_execute:
    input:
        vcf_bp=f"{config.out_dir}vcf_bp.txt",
        ref_map=f"{config.map_file}",
        bin=f"{config.root_dir}cpp/bin/interpolate-map"
    output:
        interpolated_map=f"{config.out_dir}interpolated.map"
    message:
        "Executing--interpolate map file..."
    shell:
        """
        {input.bin} {input.vcf_bp} {input.ref_map} {config.out_dir}interpolated.map;
        """

# 4.1 write segment boundary file (compile)
rule segment_boundary_file_compile:
    input:
        segment_boundary_map_cpp=f"{config.root_dir}cpp/src/segment_boundary_map.cpp",
    output:
        bin=f"{config.root_dir}cpp/bin/segment-boundary"
    message:
        "Compiling--write segment boundary file..."
    shell:
        """
        g++\
         {input.segment_boundary_map_cpp}\
         -I {config.root_dir}cpp/include/\
         -I {config.root_dir}lib/htslib/\
         -lhts\
         -o {output.bin};
        """

# 4.2 write segment boundary file (execute)
rule segment_boundary_file_execute:
    input:
        interpolated_map=f"{config.out_dir}interpolated.map",
        bin=f"{config.root_dir}cpp/bin/segment-boundary"
    output:
        segment_boundary_file=f"{config.out_dir}segment_boundary.map"
    message:
        "Executing--write segment boundary file..."
    shell:
        """
        {input.bin} {input.interpolated_map} {output.segment_boundary_file};
        """

# 5.0 slice VCF into segments
rule slice_VCF:
    input:
        vcf_file=f"{config.vcf_file}",
        segment_boundary_file=f"{config.out_dir}segment_boundary.map"
    output:
        vcf_slices_done=f"{config.out_dir}vcf_segments/vcf_slices.done"
    message:
        "Slicing VCF into segments..."
    shell:
        """
        test ! -d {config.out_dir}vcf_segments/ && mkdir {config.out_dir}vcf_segments/;
        while IFs= read -r chrm segment start_bp end_bp; do\
         echo segment ${{segment}};\
         bcftools view -h {input.vcf_file} > {config.out_dir}vcf_segments/chrm${{chrm}}.segment${{segment}}.vcf;\
         tabix {input.vcf_file} chr${{chrm}}:${{start_bp}}-${{end_bp}} >> {config.out_dir}vcf_segments/chrm${{chrm}}.segment${{segment}}.vcf;\
         bgzip {config.out_dir}vcf_segments/chrm${{chrm}}.segment${{segment}}.vcf;\
         tabix -p vcf {config.out_dir}vcf_segments/chrm${{chrm}}.segment${{segment}}.vcf.gz;\
         done < {input.segment_boundary_file};
        touch {output.vcf_slices_done};
        """
