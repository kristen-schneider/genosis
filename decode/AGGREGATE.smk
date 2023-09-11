from types import SimpleNamespace

config = SimpleNamespace(**config)

shell.prefix("""
source ~/.bashrc;
conda activate pmed;
""")

import glob
from os.path import basename

rule all:
    input:
        #f"{config.svs_results_dir}svs_results_file.txt",
        #f"{config.svs_results_dir}",
        f"{config.svs_results_dir}svs_results_file.txt",
        f"{config.svs_results_dir}",
        f"{config.cpp_bin_dir}aggregate-segments",
        f"{config.query_results_dir}segment_results.done",
        f"{config.cpp_bin_dir}aggregate-chromosomes",
        f"{config.query_results_dir}chromosome_results.done"

# 0.0 make a list of all files in sim search results dir
rule make_results_list:
    input:
        svs_results_dir=f"{config.svs_results_dir}"
    output:
        svs_results_txt=f"{config.svs_results_dir}svs_results_file.txt"
    message:
        "Writing all results file to a text file to read in..."
    shell:
        "ls {config.svs_results_dir} > {output.svs_results_txt}"
        

# 1.1 aggregate results segments (compile)
rule aggregate_segment_compile:
    input:
        file_list=f"{config.svs_results_dir}svs_results_file.txt",
        aggregate_segments_cpp=f"{config.cpp_src_dir}aggregate_segments.cpp",
        aggregate_helpers_cpp=f"{config.cpp_src_dir}aggregate_helpers.cpp",
    output:
        bin=f"{config.cpp_bin_dir}aggregate-segments"
    message:
        "Compiling--aggregating segments..."
    shell:
        "g++" \
	" {input.aggregate_segments_cpp}" \
	" {input.aggregate_helpers_cpp}" \
	" -I {config.cpp_include_dir}" \
        " -o {output.bin}"
		
# 1.2 aggregate results segments (execute)
rule aggregate_segment_execute:
    input:
        bin=f"{config.cpp_bin_dir}aggregate-segments"
    output:
        done=f"{config.query_results_dir}segment_results.done"
    message:
        "Executing--aggregating segments..."
    shell:
        "test ! -d {config.query_results_dir} && mkdir {config.query_results_dir};" \
        "echo 6. ---AGGREGATING SEGMENTS---;" \
        "{input.bin}" \
        " {config.svs_results_dir}" \
        " {config.svs_results_dir}svs_results_file.txt" \
        " {config.query_results_dir};" \
        "touch {output.done}"

# 2.1 aggregate results chromosomes
rule aggregate_chromosome_compile:
    input:
        query_results_done=f"{config.query_results_dir}segment_results.done",
        aggregate_chromosomes_cpp=f"{config.cpp_src_dir}aggregate_chromosomes.cpp",
        aggregate_helpers_cpp=f"{config.cpp_src_dir}aggregate_helpers.cpp",
    output:
        bin=f"{config.cpp_bin_dir}aggregate-chromosomes"
    message:
        "Compiling--aggregating chromosomes..."
    shell:
        "g++" \
        " {input.aggregate_chromosomes_cpp}" \
        " {input.aggregate_helpers_cpp}" \
        " -I {config.cpp_include_dir}" \
        " -o {output.bin}"

# 2.2 aggregate results chromosomes (execute)
rule aggregate_chromosomes_execute:
    input:
        bin=f"{config.cpp_bin_dir}aggregate-chromosomes"
    output:
        done=f"{config.query_results_dir}chromosome_results.done"
    message:
        "Executing--aggregating chromosomes..."
    shell:
        "echo 7. ---AGGREGATING CHROMOSOMES---;" \
        "{input.bin}" \
        " {config.query_results_dir}" \
        " {config.query_IDs};" \
        "touch {output.done}"
