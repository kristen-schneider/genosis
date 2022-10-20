from types import SimpleNamespace
configfile: "config_examples/config_GBR.yaml" # path to the config
config = SimpleNamespace(**config)

LD_LIBRARY_PATH = f"{config.conda_dir}/lib"
shell.prefix("""
set -euo pipefail;
export LD_LIBRARY_PATH=\"{LD_LIBRARY_PATH}\";
""".format(LD_LIBRARY_PATH=LD_LIBRARY_PATH))


rule all:
	input:
		f"{config.vcf_file}",
		f"{config.samples_dir}GBR_sampleIDs.txt"

# create a file with all sample IDs
# one line per sample ID
rule get_sample_IDs:
	input:
		vcf=f"{config.vcf_file}"
	output:
		sample_IDs=f"{config.samples_dir}GBR_sampleIDs.txt"
	message:
		"Creating a list of all sample IDs from VCF file..."
	shell:
		"bcftools query -l {input.vcf} > {output.sample_IDs}"
		" && touch GBR_sampleIDs.done"

