from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
#configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/AFR/AFR_config.yaml"
#configfile: "/Users/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/Users/krsc0813/chr10/config_fiji.yaml"
#configfile: "/Users/krsc0813/chr10_12/config_snakemake.yaml"
configfile: "/Users/krsc0813/AFR_pedigree/AFR_config.yaml"

config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

import glob
from os.path import basename


rule all:
    input:
        f"{config.faiss_results_dir}",
        f"{config.cpp_bin_dir}aggregate"

# 0.0 make a list of all files in sim search results dir
rule make_results_list:
    input:
        faiss_results_dir=f"{config.faiss_results_dir}"
    output:
        faiss_results_txt=f"{config.faiss_results_dir}faiss_results_file.txt"
    message:
        "Writing all results file to a text file to read in..."
    shell:
        "ls {config.faiss_results_dir} > {output.faiss_results_txt}"
        

# 1.1 aggregate results segments (compile)
rule aggregate_compile:
    input:
        main_aggregate_cpp=f"{config.cpp_src_dir}main_aggregate.cpp",
        write_query_results_cpp=f"{config.cpp_src_dir}write_query_results.cpp",
    output:
        bin=f"{config.cpp_bin_dir}aggregate"
    message:
        "Compiling--aggregating segments..."
    conda:
        f"{config.conda_pmed}"
    shell:
        "g++" \
	" {input.main_aggregate_cpp}" \
	" {input.write_query_results_cpp}" \
	" -I {config.cpp_include_dir}" \
        " -I {config.htslib_dir}" \
        " -lstdc++fs" \
        " -o {output.bin}"
		
## 1.2 encode genotypes for VCF segments (execute)
#rule encode_execute:
#    input:
#        bin=f"{config.cpp_bin_dir}encode",
#        vcf_segments=f"{config.vcf_segments_dir}{{segment}}.vcf.gz"
#    output:
#        encoding_gt=f"{config.encodings_dir}{{segment}}.gt",
#        encoding_pos=f"{config.encodings_dir}{{segment}}.pos",
#    message:
#        "Executing--encoding segments..."
#    conda:
#        f"{config.conda_pmed}"
#    shell:
#        "test ! -d {config.encodings_dir} && mkdir {config.encodings_dir};" \
#        "echo 2. ---ENCODING VCF SEGMENTS---;" \
#        "{input.bin}" \
#        " {input.vcf_segments}" \
#        " {config.root_dir}sample_IDs.txt" \
#        " {config.encoding_file}" \
#        " {config.root_dir}interpolated.map" \
#        " {config.encodings_dir};"
