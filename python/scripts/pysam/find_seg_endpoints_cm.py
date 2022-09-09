import read_map
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--map')
    parser.add_argument('--cm_max')

    return parser.parse_args()

def main():
    args = get_args()
    map_file = args.map
    max_segment_cm_length = int(args.cm_max)
    
    print('seg', 'start_bp, end_bp')
    cm_segment_endpoints_dict = read_map.find_segment_bp_endpoints_cm(map_file, max_segment_cm_length)    
    for seg in cm_segment_endpoints_dict:
        print(seg, cm_segment_endpoints_dict[seg][0], cm_segment_endpoints_dict[seg][1])

if __name__ == '__main__':
    main()
