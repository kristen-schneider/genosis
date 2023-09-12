from types import SimpleNamespace
#
config = SimpleNamespace(**config)

shell.prefix("""
#source ~/.bashrc;
. /home/sdp/miniconda3/etc/profile.d/conda.sh;
conda activate pmed;
""")

import glob
from os.path import basename

vcf_dir=f"{config.out_dir}vcf_segments/"
vcf_segments=glob.glob(vcf_dir + "*.vcf.gz")
vcf_segments=list(map(basename, vcf_segments))
vcf_segments=[".".join(v.split('.')[:-2]) for v in vcf_segments]
print(vcf_segments)
assert len(vcf_segments) > 0, "no vcf segments.."


rule all:
    input:
        expand(f"{config.out_dir}encodings/{{segment}}.gt", segment=vcf_segments),
        expand(f"{config.out_dir}encodings/{{segment}}.pos", segment=vcf_segments),

# 1.1 encode genotypes for VCF segments (compile)
rule encode_compile:
    input:
        slice_log=f"{config.out_dir}log/slice.log",
        main_encode_cpp=f"{config.root_dir}cpp/src/main_encode.cpp",
        encode_segment_cpp=f"{config.root_dir}cpp/src/encode_segment.cpp",
        read_map_cpp=f"{config.root_dir}cpp/src/read_map.cpp",
        map_encodings_cpp=f"{config.root_dir}cpp/src/map_encodings.cpp",
        utils_cpp=f"{config.root_dir}cpp/src/utils.cpp"
    output:
        bin=f"{config.root_dir}cpp/bin/encode"
    message:
        "Compiling--encoding segments..."
    shell:
        "g++" \
	" {input.main_encode_cpp}" \
	" {input.encode_segment_cpp}" \
	" {input.read_map_cpp}" \
	" {input.map_encodings_cpp}" \
	" {input.utils_cpp} " \
	" -I {config.root_dir}cpp/include/" \
        " -I {config.root_dir}lib/htslib/" \
        " -lhts" \
        " -o {output.bin}"
		
# 1.2 encode genotypes for VCF segments (execute)
rule encode_execute:
    input:
        bin=f"{config.root_dir}cpp/bin/encode",
        vcf_segments=f"{config.out_dir}vcf_segments/{{segment}}.vcf.gz"
    output:
        encoding_gt=f"{config.out_dir}encodings/{{segment}}.gt",
        encoding_pos=f"{config.out_dir}encodings/{{segment}}.pos",
    message:
        "Executing--encoding segments..."
    shell:
        "test ! -d {config.out_dir}encodings/ && mkdir {config.out_dir}encodings/;" \
        "echo 2. ---ENCODING VCF SEGMENTS---;" \
        "{input.bin}" \
        " {input.vcf_segments}" \
        " {config.out_dir}sample_IDs.txt" \
        " {config.root_dir}encoding.txt" \
        " {config.out_dir}interpolated.map" \
        " {config.out_dir}encodings/;"
