import read_map
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--map')
    parser.add_argument('--snps')

    return parser.parse_args()

def main():
    args = get_args()
    map_file = args.map
    segment_snp_length = int(args.snps)
    
    print('seg', 'start_bp, end_bp')
    snp_segment_endpoints_dict = read_map.find_segment_bp_endpoints_snp(map_file, segment_snp_length)    
    for seg in snp_segment_endpoints_dict:
        print(seg, snp_segment_endpoints_dict[seg][0], snp_segment_endpoints_dict[seg][1])

if __name__ == '__main__':
    main()
