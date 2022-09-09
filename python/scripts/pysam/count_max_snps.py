import utils
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir')
    parser.add_argument('--base_file_name')
    parser.add_argument('--chrm')
    parser.add_argument('--num_segments')
    return parser.parse_args()

def main():
    args = get_args()
    data_dir = args.data_dir
    base_file = args.data_dir + args.base_file_name
    chrm = int(args.chrm)
    num_segments = int(args.num_segments)
    
    seg_max_dict = dict()
    
    print('seg max_SNPS')
    for seg_i in range(num_segments):
        seg_file_name = base_file+str(seg_i)
        seg_max_dict = utils.read_data_file(seg_file_name)
        max_snps = max(seg_max_dict.values())
        print(seg_i, max_snps) 
    
    #for seg in range(0, len(seg_endpoint_dict)):
    #for seg in seg_endpoint_dict:
    #    start_bp = seg_endpoint_dict[seg][0]
    #    end_bp = seg_endpoint_dict[seg][1]
    #    print('evaluating segment ', seg, '. start: ', start_bp, 'end: ', end_bp)
    #    sample_snp_counts = read_vcf.sample_counts_from_vcf_block(vcf_file, chrm, start_bp, end_bp) 
    #    
    #    # open segment file and print
    #    out_file = args.out_dir + 'chr' + str(chrm) + '.cm_size.1.seg' + str(seg)
    #    o = open(out_file, 'w')
    #    o.write('sampleID SNP_count')
    #    for sample in sample_snp_counts:
    #        s = '\n' + sample + ' ' + str(sample_snp_counts[sample])
    #        o.write(s)
    #    o.close()

        #print('\n')

    #sample_snp_dict = read_vcf.sample_counts_from_vcf_block(vcf_file, chrm, 100, 61600) 

if __name__ == '__main__':
    main()
