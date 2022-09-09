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
    
    print('seg', 'start, end')
    segment_endpoints_dict = read_map.find_segment_bp_endpoints(map_file, segment_snp_length)    
    for seg in segment_endpoints_dict:
        print(seg, segment_endpoints_dict[seg])

if __name__ == '__main__':
    main()
