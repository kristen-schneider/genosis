def count_cm_per_segment(map_file=None, segment_snp_length=0):
    '''
    Returns a dictionary whose key is the segment index
    and whose value is the number of cm that each index 
    covers.

    segment length is defined by number of snps
    '''
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

def find_segment_bp_endpoints(map_file=None, segment_snp_length=0):
    '''
    Returns a dictionary whose key is the segment index
    and whose value is a tuple (start_bp, end_bp) for each index.

    segment length is defined by number of snps
    '''
    segment_bp_dict = dict()
    seg_i = 0
    snp_count = 0
    segment_bps = []

    f = open(map_file, 'r')
    for line in f:
        L = line.strip().split()
        snp_count += 1
        chrm = L[0]
        cm = float(L[2])
        bp = int(L[3])
        segment_bps.append(bp)
        
        if snp_count == segment_snp_length:
            segment_bp_dict[seg_i] = (segment_bps[0], segment_bps[-1])
            seg_i += 1
            snp_count = 0
    f.close()

    return segment_bp_dict


def count_snps_per_cm(map_file=None, cm_max=0):
    '''
    Returns a dictionary whose key is the segment index
    and whose value is the number of SNPs that each index 
    covers.

    segment length is defined by a fixed cm length
    '''

    cm_snps_dict = dict()
    seg_cm_start = 0
    seg_i = 0
    snp_count = 0
    
    max_cm_length = cm_max

    f = open(map_file, 'r')
    for line in f:
        L = line.strip().split()
        snp_count += 1
        chrm = L[0]
        cm = float(L[2])
        bp = int(L[3])
        
        
        snp_count += 1
        seg_cm_length = cm - seg_cm_start
        if  seg_cm_length >= max_cm_length:
            seg_cm_start = cm
            cm_snps_dict[seg_i] = snp_count
            seg_i += 1

    seg_i += 1
    cm_snps_dict[seg_i] = snp_count
    f.close()

    return cm_snps_dict


def make_bp_cm_dict(map_file=None):
    '''
    Returns a dictionary whose key is the 
    basepair start pos and whose value 
    is the cm distance at that basepair.
    '''
    bp_cm_dict = dict()
    
    f = open(map_file, 'r')
    for line in f:
        L = line.strip().split()
        chrm = L[0]
        cm = float(L[2])
        bp = int(L[3])

        bp_cm_dict[bp] = cm

    return bp_cm_dict
