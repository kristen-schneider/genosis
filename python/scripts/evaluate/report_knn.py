import argparse
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query_IDs')
    parser.add_argument('--ss_sample_results_dir')
    parser.add_argument('--out_dir')

    return parser.parse_args()

def main():
    args = get_args()
    query_IDs_txt = args.query_IDs
    ss_sample_results_dir = args.ss_sample_results_dir + '/'
    out_dir = args.out_dir + '/'

    query_IDs = read_query_IDs(query_IDs_txt)
    for q_ID in query_IDs:
        
        queryID_hap0 = ss_sample_results_dir + q_ID + '_0/' + 'all_chromosomes.csv'
        hap0_results = read_hap_results(queryID_hap0)

        queryID_hap1 = ss_sample_results_dir + q_ID + '_1/' + 'all_chromosomes.csv'
        hap1_results = read_hap_results(queryID_hap1)

        write_combined_haps(q_ID, hap0_results, hap1_results, out_dir)

def write_combined_haps(q_ID, hap0_results, hap1_results, out_dir):

    print(q_ID, hap0_results, hap1_results)

    svs_csv = out_dir + q_ID + '.svs.csv'
    svs_file = open(svs_csv, 'w')
    
    pop_csv = out_dir + q_ID + '.pop.csv'
    pop_file = open(pop_csv, 'w')
    
    ibd_csv = out_dir + q_ID + '.ibd.csv'
    ibd_file = open(ibd_csv, 'w')

    svs_file.write(q_ID + '_0\n')
    pop_file.write(q_ID + '_0\n')
    ibd_file.write(q_ID + '_0\n')
    for match_id in hap0_results[q_ID]:
        svs_file.write(' ' + match_id + ': ' + str(hap0_results[q_ID][match_id][0]) + '\n')
        pop_file.write(' ' + match_id + ': ' + str(hap0_results[q_ID][match_id][1]) + '\n')
        ibd_file.write(' ' + match_id + ': ' + str(hap0_results[q_ID][match_id][2]) + '\n')

    svs_file.write(q_ID + '_1\n')
    pop_file.write(q_ID + '_1\n')
    ibd_file.write(q_ID + '_1\n')
    for match_id in hap1_results[q_ID]:
        svs_file.write(' ' + match_id + ': ' + str(hap1_results[q_ID][match_id][0]) + '\n')
        pop_file.write(' ' + match_id + ': ' + str(hap1_results[q_ID][match_id][1]) + '\n')
        ibd_file.write(' ' + match_id + ': ' + str(hap1_results[q_ID][match_id][2]) + '\n')


def read_hap_results(hap_results_csv):
    h_file = open(hap_results_csv, 'r')
    hap_results = defaultdict(dict)
    sample_hap = h_file.readline().strip()
    sample = sample_hap[:-2]
    header = h_file.readline()

    for line in h_file:
        line = line.strip().split(',')
        chrm = int(line[0])
        match_hap = line[1]
        match_id = match_hap
        svs_score = float(line[2])
        pop_score = float(line[3])
        merge_gap_score = float(line[6])
        try:
            hap_results[sample][match_id][0] += svs_score
            hap_results[sample][match_id][1] += pop_score
            hap_results[sample][match_id][2] += merge_gap_score
        except KeyError:
            hap_results[sample][match_id] = [svs_score, pop_score, merge_gap_score]

    return hap_results


def read_query_IDs(query_IDs_txt):
    query_IDs = open(query_IDs_txt).read().splitlines()
    return query_IDs

if __name__ == '__main__':
    main()
