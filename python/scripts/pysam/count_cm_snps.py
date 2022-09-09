import argparse
import read_map

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--map')
    parser.add_argument('--cm_max')

    return parser.parse_args()

def main():
    args = get_args()
    map_file = args.map
    max_segment_cm_length = int(args.cm_max)
    
    print('seg', 'SNPs')
    cm_snps_dict = read_map.count_snps_per_cm(map_file, max_segment_cm_length)
    for seg in cm_snps_dict:
        print(seg, cm_snps_dict[seg])

    


if __name__ == '__main__':
    main()
