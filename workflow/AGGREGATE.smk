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

import glob
from os.path import basename

rule all:
    input:
        svs_results_txt=f"{config.out_dir}svs_results/svs_results_file.txt",
        top_hits=f"{config.out_dir}TOP_HITS.txt"
        #f"{config.out_dir}svs_sample_results/segment_results.done",
        #f"{config.out_dir}svs_sample_results/chromosome_results.done"
        #f"{config.out_dir}svs_sample_results/knn_summary.done"

# 0.0 make a list of all files in sim search results dir
rule make_results_list:
    input:
        svs_results_dir=f"{config.out_dir}svs_results/"
    output:
        svs_results_txt=f"{config.out_dir}svs_results/svs_results_file.txt"
    message:
        "Writing all results file to a text file to read in..."
    shell:
        """
        ls {config.out_dir}svs_results/ > {output.svs_results_txt};
        """

# 1.1 aggregate svs_results with ryan's agg scripts
rule aggregate_ryan:
    input:
        ryan_agg=f"{config.root_dir}python/scripts/aggregate/ryan_agg.py",
    output:
        top_hits=f"{config.out_dir}TOP_HITS.txt"
    message:
        "Running Ryan's aggregation script..."
    shell:
        """
        python {config.root_dir}python/scripts/aggregate/ryan_agg.py \
          --id_file {config.out_dir}query_hap_IDs.txt \
          --data_dir {config.out_dir}svs_results/ \
          --num_threads 20\
          --top_N 20\
          --out_file {output.top_hits}
        """

## 1.1 aggregate results segments (compile)
#rule aggregate_segment_compile:
#    input:
#        file_list=f"{config.out_dir}svs_results/svs_results_file.txt",
#        aggregate_segments_cpp=f"{config.root_dir}cpp/src/aggregate_segments.cpp",
#        aggregate_helpers_cpp=f"{config.root_dir}cpp/src/aggregate_helpers.cpp",
#    output:
#        bin=f"{config.root_dir}cpp/bin/aggregate-segments"
#    message:
#        "Compiling--aggregating segments..."
#    shell:
#        """
#        g++\
#	 {input.aggregate_segments_cpp}\
#	 {input.aggregate_helpers_cpp}\
#	 -I {config.root_dir}cpp/include/\
#         -o {output.bin};
#        """
#		
## 1.2 aggregate results segments (execute)
#rule aggregate_segment_execute:
#    input:
#        bin=f"{config.root_dir}cpp/bin/aggregate-segments"
#    output:
#        done=f"{config.out_dir}svs_sample_results/segment_results.done"
#    message:
#        "Executing--aggregating segments..."
#    shell:
#        """
#        test ! -d {config.out_dir}svs_sample_results/ && mkdir {config.out_dir}svs_sample_results/;
#        {input.bin}\
#         {config.out_dir}svs_results/\
#         {config.out_dir}svs_results/svs_results_file.txt\
#         {config.out_dir}svs_sample_results/
#        """
#
## 2.1 aggregate results chromosomes (compile)
#rule aggregate_chromosome_compile:
#    input:
#        query_results_done=f"{config.out_dir}svs_sample_results/segment_results.done",
#        aggregate_chromosomes_cpp=f"{config.root_dir}cpp/src/aggregate_chromosomes.cpp",
#        aggregate_helpers_cpp=f"{config.root_dir}cpp/src/aggregate_helpers.cpp",
#    output:
#        bin=f"{config.root_dir}cpp/bin/aggregate-chromosomes"
#    message:
#        "Compiling--aggregating chromosomes..."
#    shell:
#        """
#        g++\
#         {input.aggregate_chromosomes_cpp}\
#         {input.aggregate_helpers_cpp}\
#         -I {config.root_dir}cpp/include/\
#         -o {output.bin};
#        """
#
## 2.2 aggregate results chromosomes (execute)
#rule aggregate_chromosomes_execute:
#    input:
#        bin=f"{config.root_dir}cpp/bin/aggregate-chromosomes"
#    output:
#        done=f"{config.out_dir}svs_sample_results/chromosome_results.done"
#    message:
#        "Executing--aggregating chromosomes..."
#    shell:
#        """
#        {input.bin}\
#         {config.out_dir}svs_sample_results/\
#         {config.out_dir}query_IDs.txt;
#        """
#
### 3.0 report knn for samples
##rule report sample_knn:
##    input:
##        query_IDs=f"{config.out_dir}query_IDs.txt",
##        chrm=f"{config.out_dir}svs_sample_results/chromosome_results.done"
##    output:
##        knn_summary=f"{config.out_dir}svs_sample_results/knn_summary.done"
##    message:
##        "Writing a summary result file with knn for all samples..."
##    shell:
##        """
##        python {config.root_dir}python/scripts/evaluate/report_knn.py\
##         --query_IDs {input.query_IDs}\
##         --ss_sample_results_dir {config.out_dir}svs_sample_results/\
##         --out_dir {config.out_dir}svs_sample_results/;
##         touch {output.knn_summary};
##        """
