import read_vcf
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--chrm')
    parser.add_argument('--cm_endpoints')
    parser.add_argument('--out_dir')
    return parser.parse_args()

def main():
    args = get_args()
    vcf_file = args.vcf
    chrm = args.chrm
    cm_endpoint_file = args.cm_endpoints
    
    seg_freq_dict = dict()
    
    seg_endpoint_dict = read_vcf.read_endpoint_file(cm_endpoint_file)
    for seg in range(0, len(seg_endpoint_dict)):
    #for seg in seg_endpoint_dict:
        start_bp = seg_endpoint_dict[seg][0]
        end_bp = seg_endpoint_dict[seg][1]
        print('evaluating segment ', seg, '. start: ', start_bp, 'end: ', end_bp)
        sample_snp_counts = read_vcf.sample_counts_from_vcf_block(vcf_file, chrm, start_bp, end_bp) 
        
        # open segment file and print
        out_file = args.out_dir + 'chr' + str(chrm) + '.cm_size.1.seg' + str(seg)
        o = open(out_file, 'w')
        o.write('sampleID SNP_count')
        for sample in sample_snp_counts:
            s = '\n' + sample + ' ' + str(sample_snp_counts[sample])
            o.write(s)
        o.close()

        #print('\n')

    #sample_snp_dict = read_vcf.sample_counts_from_vcf_block(vcf_file, chrm, 100, 61600) 

if __name__ == '__main__':
    main()
