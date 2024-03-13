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

# get a list of all vcf segments
VCF_DIR=f"{config.out_dir}vcf_segments/"
VCF_SEGMENTS=glob.glob(VCF_DIR + "*.vcf.gz")
VCF_SEGMENTS=list(map(basename, VCF_SEGMENTS))
VCF_SEGMENTS=[vcf_seg.replace('.vcf.gz', '') for vcf_seg in VCF_SEGMENTS]
assert len(VCF_SEGMENTS) > 0, "no vcf segments.."

rule all:
    input:
        gt_encodings=expand(f"{config.out_dir}encodings/{{segment}}.gt", segment=VCF_SEGMENTS),
	zeros=expand(f"{config.out_dir}encodings/{{segment}}.done", segment=VCF_SEGMENTS),
        database_hap_IDs=f"{config.out_dir}database_hap_IDs.txt",
        query_hap_IDs=f"{config.out_dir}query_hap_IDs.txt"
 
# 1.1 encode genotypes for VCF segments (compile)
rule encode_compile:
    input:
        main_encode_cpp=f"{config.root_dir}cpp/src/main_encode.cpp",
        encode_segment_cpp=f"{config.root_dir}cpp/src/encode_segment.cpp",
        read_map_cpp=f"{config.root_dir}cpp/src/read_map.cpp",
        map_encodings_cpp=f"{config.root_dir}cpp/src/map_encodings.cpp",
        utils_cpp=f"{config.root_dir}cpp/src/utils.cpp"
    output:
        encode_bin=f"{config.root_dir}cpp/bin/encode"
    message:
        "Compiling--encoding segments..."
    shell:
        """
        g++\
	 {input.main_encode_cpp}\
	 {input.encode_segment_cpp}\
	 {input.read_map_cpp}\
	 {input.map_encodings_cpp}\
	 {input.utils_cpp}\
	 -I {config.root_dir}cpp/include/\
         -I {config.root_dir}lib/htslib/\
         -lhts\
         -o {output.encode_bin}
        """
	
# 1.2 encode genotypes for VCF segments (execute)
rule encode_execute:
    input:
        encode_bin=f"{config.root_dir}cpp/bin/encode",
        vcf_segments=f"{config.out_dir}vcf_segments/{{segment}}.vcf.gz"
    output:
        gt_encodings=f"{config.out_dir}encodings/{{segment}}.gt",
        pos_encodings=f"{config.out_dir}encodings/{{segment}}.pos"
    message:
        "Executing--encoding segments..."
    shell:
        """
        test ! -d {config.out_dir}encodings/ && mkdir {config.out_dir}encodings/;\
        {input.encode_bin}\
         {input.vcf_segments}\
         {config.out_dir}sample_IDs.txt\
         {config.root_dir}encoding.txt\
         {config.out_dir}interpolated.map\
         {config.out_dir}encodings/;
        """

# 2.0 remove empty entries
rule remove_empty_encodings:
    input:
        pos_encodings=f"{config.out_dir}encodings/{{segment}}.pos"
    output:
        zeros=f"{config.out_dir}encodings/{{segment}}.done"
    message:
        "Removing positional encodings with empty sample entries"
    shell:
        """
        python {config.root_dir}python/scripts/check_pos_encodings.py\
         --pos_enc {input.pos_encodings} > {output.zeros}
        """

# 3.0 get haplotype IDs for database samples
rule hap_IDs:
    input:
        hap_script=f"{config.root_dir}helper_scripts/make_hap_IDs.sh",
        sample_IDs=f"{config.out_dir}sample_IDs.txt",
        database_IDs=f"{config.database_IDs}",
        query_IDs=f"{config.query_IDs}"
    output:
        sample_hap_IDs=f"{config.out_dir}sample_hap_IDs.txt",
        database_hap_IDs=f"{config.out_dir}database_hap_IDs.txt",
        query_hap_IDs=f"{config.out_dir}query_hap_IDs.txt"
    message:
        "Getting a list of all sampleIDs, databaseIDs, and queryIDs"
    shell:
        """
        bash {input.hap_script} {input.sample_IDs} > {output.sample_hap_IDs};
        bash {input.hap_script} {input.database_IDs} > {output.database_hap_IDs};
        bash {input.hap_script} {input.query_IDs} > {output.query_hap_IDs};
        """
