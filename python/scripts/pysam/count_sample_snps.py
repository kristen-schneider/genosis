import read_vcf
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--chrm')
    parser.add_argument('--cm_endpoints')

    return parser.parse_args()

def main():
    args = get_args()
    vcf_file = args.vcf
    chrm = args.chrm
    cm_endpoint_file = args.cm_endpoints
    
    #seg_endpoint_dict = read_vcf.read_endpoint_file(cm_endpoint_file)
    #for seg in seg_endpoint_dict:
    #    start_bp = seg_endpoint_dict[seg][0]
    #    end_bp = seg_endpoint_dict[seg][1]
    #    sample_snp_counts = read_vcf.sample_counts_from_vcf_block(vcf_file, chrm, start_bp, end_bp) 

    read_vcf.sample_counts_from_vcf_block(vcf_file, chrm, 100, 61600) 

if __name__ == '__main__':
    main()
