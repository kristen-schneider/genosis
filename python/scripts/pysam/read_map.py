def count_cm_per_segment(map_file=None, segment_snp_length=0):
    
    segment_cm_dict = dict()
    seg_i = 0
    snp_count = 0
    start_cm = 0

    f = open(map_file, 'r')
    for line in f:
        L = line.strip().split()
        snp_count += 1
        chrm = L[0]
        cm = float(L[2])
        bp = int(L[3])

        if snp_count == segment_snp_length:
            cm_dist = cm - start_cm
            segment_cm_dict[seg_i] = cm_dist
            seg_i += 1
            start_cm = cm
            snp_count = 0 
    f.close()

    return segment_cm_dict

def make_bp_cm_dict(map_file=None):
    bp_cm_dict = dict()
    
    f = open(map_file, 'r')
    for line in f:
        L = line.strip().split()
        chrm = L[0]
        cm = float(L[2])
        bp = int(L[3])

        bp_cm_dict[bp] = cm

    return bp_cm_dict
