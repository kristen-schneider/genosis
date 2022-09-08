import pysam

def read_vcf_block(vcf_file, start_bp, end_bp):
    
    vcf_f = pysam.VariantFile(vcf_file)
    
    for snp in vcf_f.fetch(start_bp, end_bp):
        print(snp)
