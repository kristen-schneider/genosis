import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--map')
    parser.add_argument('--out_dir')
    parser.add_argument('--num_snps')

    return parser.parse_args()

def main():
    args = get_args()
    bp_start_end = dict()
    cm_map_dict = dict()
    cm_est_count = dict()
    cm_map_count = dict()

    bp_start_end = get_all_bp_count(args.vcf, int(args.num_snps))
    cm_map_dict = get_cm_map(args.map)
    cm_est_count = get_segment_cm_estimated(bp_start_end)
    cm_map_count = get_segment_cm_map(bp_start_end, cm_map_dict)

    bp_distance_file = args.out_data + str(args.num_snps) + '.basepair_distance.txt'
    write_snp_data(bp_start_end, bp_distance_file)
    cm_est_distance_file = args.out_data + str(args.num_snps) + '.centimorgan_est_distance.txt'
    write_snp_data(cm_est_count, cm_est_distance_file)
    cm_map_distance_file = args.out_data + str(args.num_snps) + '.centimorgan_map_distance.txt'
    write_snp_data(cm_map_count, cm_map_distance_file)

def write_snp_data(snp_data, out_file):
    f = open(out_file, 'w')
    for seg in snp_data:
        dist = snp_data[seg]
        if isinstance(dist, list):
            dist = dist[1] - dist[0]
        line = str(seg) + ',' + str(dist) + '\n'
        f.write(line)
    f.close()

def get_all_bp_count(vcf_file, num_snps):
    bp_start_end = dict()
    seg_i = 0
    seg_snp_count = 0
    seg_positions = []
    f = open(vcf_file)
    for line in f:
        if '#' in line:
            continue
        elif (seg_snp_count < num_snps):
            pos = int(line.strip().split()[1])
            seg_positions.append(pos)
            seg_snp_count += 1
        else:
            # save off last segment
            bp_start_end[seg_i] = [seg_positions[0], seg_positions[-1]]
            seg_i += 1
            # start new segment
            pos = int(line.strip().split()[1])
            seg_positions = [pos]
            seg_snp_count = 1
    f.close()
    # anything leftover
    if (seg_snp_count <= num_snps):
        seg_i += 1
        bp_start_end[seg_i] = [seg_positions[0], seg_positions[-1]]

    return bp_start_end

def get_cm_map(interpolated_map):
    pos_cm_dict = dict()
    f = open(interpolated_map, 'r')
    for line in f:
        L = line.strip().split()
        cm = L[2]
        pos = L[3]
        pos_cm_dict[int(pos)] = float(cm)
    f.close()
    return pos_cm_dict

def get_segment_cm_estimated(bp_start_end):
    # 1e6 bp = 1cM
    cm_est_count = dict()
    for seg in bp_start_end:
        bp_distance = bp_start_end[seg][1] - bp_start_end[seg][0]
        bp_per_cm = 1e6
        cm_estimate = bp_distance / bp_per_cm
        cm_est_count[seg] = cm_estimate
    return cm_est_count

def get_segment_cm_map(bp_start_end, cm_map_dict):
    # cM distance from map file
    cm_map_count = dict()
    for seg in bp_start_end:
        startPOS = bp_start_end[seg][0]
        endPOS = bp_start_end[seg][1]
        startCM = cm_map_dict[startPOS]
        endCM = cm_map_dict[endPOS]
        cm_map = endCM-startCM
        cm_map_count[seg] = cm_map
    return cm_map_count

if __name__ == '__main__':
    main()
