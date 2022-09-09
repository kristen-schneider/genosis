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
    
    print('seg', 'cM')
    segment_cm_dict = read_map.count_cm_per_segment(map_file, segment_snp_length)    
    for seg in segment_cm_dict:
        print(seg, segment_cm_dict[seg])

if __name__ == '__main__':
    main()
