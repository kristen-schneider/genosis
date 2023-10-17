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

# 1.0 make relations file (if pedigree)
rule label_relations:
    input:
        ped_file=f"{config.ped}",
        root_file=f"{config.roots}",
    message:
        "Labeling relationships for pedigree data..."
    output:
        relations_file=f"{config.out_dir}samples.relations"
    shell:
        "python {config.root_dir}python/scripts/evaluate/relations.py" \
        " --ped {input.ped_file}" \
        " --root_file {input.root_file}" \
        " --relations {output.relations_file}" \
        " --gen 6;"

# plot KNN results
# read ground truth IBD
# plot KNN vs ground truth


# 2.0 plot summmary data
rule plot_summary:
    input:
        knn_summary=f"{config.out_dir}svs_sample_results/knn_summary.done",
	relations=f"{config.out_dir}samples.relations",
	ancestry=f"{config.ancestry}",
        query_IDs=f"{config.out_dir}query_IDs.txt",
        svs_results_dir=f"{config.out_dir}svs_sample_results/"
    output:
        plot_summary=f"{config.out_dir}svs_sample_results/plot_summary.done"
    message:
        "Plotting summary results..."
    shell:
        "python {config.root_dir}python/scripts/evaluate/label_category.py" \
	" --relations {input.relations}" \
        " --ancestry {input.ancestry}" \
	" --samples {input.sample_IDs}" \
	" --data_dir {input.svs_results_dir} > {output.plot_summary};"
