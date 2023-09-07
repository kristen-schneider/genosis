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
		f"{config.log_dir}plink2_kingship.log",
		f"{config.log_dir}phased_ibd_vcf.log",
		f"{config.log_dir}ibd_map.log",
		f"{config.log_dir}phased_ibd.log",
		f"{config.log_dir}ped.log"

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

# 3.1 bgzip decompress vcf file for phased ibd
rule phased_ibd_vcf:
	input:
		vcf_file=f"{config.vcf_file}"
	log:
		pibd_vcf_log=f"{config.log_dir}phased_ibd_vcf.log"
	message:
		"bgzip decompress vcf file for phased ibd..."
	conda:
		f"{config.conda_pmed}"
	shell:
		"bgzip -d -k {input.vcf_file};"

# 3.2 rewrite map file for phased ibd and ilash
rule phased_ibd_map:
	input:
		map_file=f"{config.root_dir}interpolated.map"
	output:
		ibd_map=f"{config.root_dir}ibd.map"
	log:
		pibd_map_log=f"{config.log_dir}ibd_map.log"
	message:
		"making new map file in format for phased ibd and ilash..."
	conda:
		f"{config.conda_pmed}"
	shell:
		"awk '{{print $1, $2, $2, $3}}' {input.map_file} > {output.ibd_map}"	

# 3.3 run phased ibd
rule phased_ibd:
	input:
		vcf_file=f"{config.vcf_file}",
		map_file=f"{config.root_dir}ibd.map"
	log:
		phased_ibd=f"{config.log_dir}phased_ibd.log"
	benchmark:
		f"{config.benchmark_dir}phased_ibd.tsv"
	message:
		"Running phased IBD..."
	conda:
		f"{config.conda_pmed}"
	shell:
		"filename=$(basename {input.vcf_file});" \
		"decompressed_vcf=${{filename%.*}};" \
		"python {config.python_dir}ibd/phased_ibd.py" \
		"	--vcf $decompressed_vcf" \
		"	--map {input.map_file}" \
		"	--out {config.compare_dir}phased_ibd.csv > {log.phased_ibd}"
		 
# 4.1 ped file for iLASH 
rule ped_file:
	input:
		vcf_file=f"{config.vcf_file}"
	log:
		plink_ped=f"{config.log_dir}ped.log"
	benchmark:
		f"{config.benchmark_dir}ped.tsv"
	message:
		"Creating ped file..."
	conda:
		f"{config.conda_pmed}"
	shell:
		"plink --vcf {input.vcf_file}" \
		"	--recode 01" \
		"	--output-missing-genotype ." \
		"	--out {config.compare_dir}ex-ped"

# 4.2 Run iLASH 
rule run_ilash:
	input:
		vcf_file=f"{config.vcf_file}",
		ilash_config=f"{config.ilash_config}"
	log:
		ilash_log=f"{config.log_dir}ilash.log"
	benchmark:
		f"{config.benchmark_dir}ilash.tsv"
	message:
		"Running iLASH..."
	conda:
		f"{config.conda_pmed}"
	shell:
		"../iLASH/build/ilash {ilash.config}" \
