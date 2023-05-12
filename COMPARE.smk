from types import SimpleNamespace
#configfile: "/home/sdp/pmed-local/data/1KG/config_snakemake.yaml"
configfile: "/home/sdp/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/precision-medicine/example/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/1kg/config_snakemake.yaml"
#configfile: "/scratch/alpine/krsc0813/data/SAS/SAS_config.yaml"
config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_pmed}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))

rule all:
	input:
		f"{config.log_dir}plink_genome.log",
		f"{config.log_dir}plink2_kingship.log"

# 1.0 plink --genome
rule plink_genome:
	input:
		vcf_file=f"{config.vcf_file}"
	log:
		plink_genome_log=f"{config.log_dir}plink_genome.log"
	benchmark:
		f"{config.benchmark_dir}plink_genome.tsv"
	message:
		"Running plink genome..."
	conda:
		f"{config.conda_pmed}"
	shell:
		"plink --vcf {input.vcf_file} " \
		"	--genome " \
		"	--out {config.compare_dir}ex-plink"

# 2.0 plink2 --kingship
rule plink2_kingship:
	input:
		vcf_file=f"{config.vcf_file}"
	log:
		plink2_kingship_log=f"{config.log_dir}plink2_kingship.log"
	benchmark:
		f"{config.benchmark_dir}plink2_kingship.tsv"
	message:
		"Running plink2 kingship..."
	conda:
		f"{config.conda_pmed}"
	shell:
		"plink2 --vcf {input.vcf_file} " \
		"	--make-king-table " \
		"	--out {config.compare_dir}ex-plink2"

