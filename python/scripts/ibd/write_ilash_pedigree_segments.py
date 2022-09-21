import argparse
import os

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ilash_file')
    parser.add_argument('--query')
    parser.add_argument('--dad')
    parser.add_argument('--mom')
    parser.add_argument('--hap')
    parser.add_argument('--map')
    parser.add_argument('--cm')
    return parser.parse_args()

def main():
    args = get_args()
    # get start and end bp of each segment
    map_cm_segments = make_map_cm_segments(args.map, int(args.cm))
    #print(len(map_cm_segments), map_cm_segments)

    # get dictionary of matches
    query_hap_dict = read_ilash_segments(args.ilash_file, args.query, args.dad, args.mom, args.hap)
    #for m in query_hap_dict:
    #    print(m, query_hap_dict[m])
    
    # get cm segment for each match
    match_to_cm_segment(query_hap_dict, map_cm_segments)

def match_to_cm_segment(query_hap_dict, map_cm_segments):
    print('sampleID, start_bp, end_bp, start_seg, end_seg')
    for match in query_hap_dict:
        for start_end_bp in query_hap_dict[match]:
            start = start_end_bp[0]
            end = start_end_bp[1]
            start_end_seg = linear_search(map_cm_segments, start, end)
            print(match, start_end_bp[0], start_end_bp[1], start_end_seg[0], start_end_seg[1])

def linear_search(map_cm_segments, start, end):
    i = 0
    for start_end in map_cm_segments:
        if start >= start_end[0] and start <= start_end[1]:
            start_i = i
        if end <= start_end[1]:
            end_i = i
            return(start_i, end_i)
        i += 1
    end_i = i
    return(start_i, end_i)

def make_map_cm_segments(map_file, cm_slice):
    cm_start_end_bp = []

    mf = open(map_file, 'r')
    cm_max = cm_slice
    first_line = mf.readline()
    L = first_line.strip().split()
    if float(L[2]) <= cm_slice:
        start_bp = int(L[3])

    for line in mf:
        L = line.strip().split()
        cm = float(L[2])
        bp = int(L[3])
        
        if cm >= cm_max:
            cm_max = cm + cm_slice
            end_bp = int(L[3])
            segment = (start_bp, end_bp)
            cm_start_end_bp.append(segment)
            start_bp = end_bp + 1
    
    return cm_start_end_bp
        
def read_ilash_segments(ilash_file, query, dad, mom, hap):
    query_hap_dict = dict()
    query_hap = query+'_'+hap

    f = open(ilash_file, 'r')
    for line in f:
        L = line.strip().split()
        sample1 = L[0]
        sample1_hap = L[1]
        sample2 = L[2]
        sample2_hap = L[3]
        bp_start = int(L[5])
        bp_end = int(L[6])
        
        if query_hap == sample1_hap:
            if dad == sample2 or mom == sample2 or query == sample2:
                try:
                    query_hap_dict[sample2_hap].append((bp_start, bp_end))
                except KeyError:
                    query_hap_dict[sample2_hap] = [(bp_start, bp_end)]
        elif query_hap == sample2_hap:
            if dad == sample1 or mom == sample1 or query == sample1:
                try:
                    query_hap_dict[sample1_hap].append((bp_start, bp_end))
                except KeyError:
                    query_hap_dict[sample1_hap] = [(bp_start, bp_end)]
    f.close()
    return query_hap_dict

if __name__ == '__main__':
    main()
