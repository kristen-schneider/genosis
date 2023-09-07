from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/AFR/AFR_config.yaml"
configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/Users/krsc0813/chr10/config_fiji.yaml"
#configfile: "/Users/krsc0813/chr15_20/config_snakemake.yaml"
#configfile: "/Users/krsc0813/AFR_pedigree/config_AFR.yaml"

config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

import glob
from os.path import basename

vcf_dir=f"{config.vcf_segments_dir}"
vcf_segments=glob.glob(vcf_dir + "*.vcf.gz")
vcf_segments=list(map(basename, vcf_segments))
vcf_segments=[".".join(v.split('.')[:-2]) for v in vcf_segments]
print(vcf_segments)
assert len(vcf_segments) > 0, "no vcf segments.."


rule all:
    input:
        expand(f"{config.encodings_dir}{{segment}}.gt", segment=vcf_segments),
        expand(f"{config.encodings_dir}{{segment}}.pos", segment=vcf_segments),

# 1.1 encode genotypes for VCF segments (compile)
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
		
# 1.2 encode genotypes for VCF segments (execute)
rule encode_execute:
    input:
        bin=f"{config.cpp_bin_dir}encode",
        vcf_segments=f"{config.vcf_segments_dir}{{segment}}.vcf.gz"
    output:
        encoding_gt=f"{config.encodings_dir}{{segment}}.gt",
        encoding_pos=f"{config.encodings_dir}{{segment}}.pos",
    message:
        "Executing--encoding segments..."
    conda:
        f"{config.conda_pmed}"
    shell:
        "test ! -d {config.encodings_dir} && mkdir {config.encodings_dir};" \
        "echo 2. ---ENCODING VCF SEGMENTS---;" \
        "{input.bin}" \
        " {input.vcf_segments}" \
        " {config.root_dir}sample_IDs.txt" \
        " {config.encoding_file}" \
        " {config.root_dir}interpolated.map" \
        " {config.encodings_dir};"
