import argparse
from collections import defaultdict

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_IDs')
    parser.add_argument('--ss_sample_results_dir')
    parser.add_argument('--out_dir')

    return parser.parse_args()

def main():
    args = get_args()
    sample_IDs_txt = args.sample_IDs
    ss_sample_results_dir = args.ss_sample_results_dir
    out_dir = args.out_dir

    sample_IDs = read_sample_IDs(sample_IDs_txt)
    for sample_ID in sample_IDs:
        print(sample_ID)
        sampleID_hap0 = ss_sample_results_dir + sample_ID + '_0/' + 'all_chromosomes.csv'
        hap0_results = read_hap_results(sampleID_hap0)

        sampleID_hap1 = ss_sample_results_dir + sample_ID + '_1/' + 'all_chromosomes.csv'
        hap1_results = read_hap_results(sampleID_hap1)

        # combined_haps = combine_key_haps(hap0_results, hap1_results)

        write_combined_haps(sample_ID, hap0_results, hap1_results, out_dir)

def write_combined_haps(sample_ID, hap0_results, hap1_results, out_dir):

    out_csv = out_dir + sample_ID + '.csv'
    out_file = open(out_csv, 'w')
    out_file.write(sample_ID + '_0\n')
    for match_id in hap0_results[sample_ID]:
        out_file.write(' ' + match_id + ': ' + str(hap0_results[sample_ID][match_id]) + '\n')

    out_file.write(sample_ID + '_1\n')
    for match_id in hap1_results[sample_ID]:
        out_file.write(' ' + match_id + ': ' + str(hap1_results[sample_ID][match_id]) + '\n')


def combine_key_haps(hap0_results, hap1_results):
    combined_haps = defaultdict(dict)
    for sample in hap0_results:
        for match_id in hap0_results[sample]:
            combined_haps[sample][match_id] = (
                    hap0_results[sample][match_id] + hap1_results[sample][match_id])
    return combined_haps

def read_hap_results(hap_results_csv):
    h_file = open(hap_results_csv, 'r')
    hap_results = defaultdict(dict)
    sample_hap = h_file.readline().strip()
    sample = sample_hap.split('_')[0]
    header = h_file.readline()

    for line in h_file:
        line = line.strip().split(',')
        chrm = int(line[0])
        match_hap = line[1]
        match_id = match_hap.split('_')[0]
        ibd_score = float(line[5])
        try:
            hap_results[sample][match_id] += ibd_score
        except KeyError:
            hap_results[sample][match_id] = ibd_score

    return hap_results


def read_sample_IDs(sample_IDs_txt):
    sample_IDs = open(sample_IDs_txt).read().splitlines()
    return sample_IDs

if __name__ == '__main__':
    main()
