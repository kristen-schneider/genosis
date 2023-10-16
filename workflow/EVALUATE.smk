from types import SimpleNamespace
#
config = SimpleNamespace(**config)

#shell.prefix("""
#. /opt/conda/etc/profile.d/conda.sh
#conda activate pmed;
#""")
shell.prefix("""
source  ~/.bashrc
conda activate pmed;
""")



import glob
from os.path import basename

rule all:
    input:
        f"{config.out_dir}svs_sample_results/plot_summary.done"

# make relations file
# plot KNN results
# read ground truth IBD
# plot KNN vs ground truth


# 1.0 plot summmary data
rule plot_summary:
    input:
        knn_summary=f"{config.out_dir}svs_sample_results/knn_summary.done",
	relations=f"{config.out_dir}samples.relations",
        query_IDs=f"{config.out_dir}query_IDs.txt",
        svs_results_dir=f"{config.out_dir}svs_sample_results/"
    output:
        plot_summary=f"{config.out_dir}svs_sample_results/plot_summary.done"
    message:
        "Plotting summary results..."
    shell:
        "python {config.root_dir}python/scripts/decode/lump_relations.py" \
	" --relations {input.relations}" \
	" --samples {input.sample_IDs}" \
	" --data_dir {input.svs_results_dir};" \
	" touch {output.plot_summary};" 
