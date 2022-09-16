import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pedigree_data')
    parser.add_argument('--samples_dir')
    return parser.parse_args()

def main():
    args = get_args()
    make_fam_summary_file(args.pedigree_data, args.samples_dir)
    #pop_samples = make_pop_samples_dict(args.sample_data)

def make_fam_summary_file(pedigree_file, out_dir):
    o = open(out_dir + 'pedigree.summary', 'w')

    f = open(pedigree_file, 'r')
    header = f.readline().split()
    o.write(header[0] + ' ' + header[1] + ' ' + header[2] + '\n')
    for line in f:
        L = line.strip().split()
        sample_ID = L[0]
        father_ID = L[1]
        mother_ID = L[2]
        
        if father_ID != '0' or mother_ID != '0':
            out_line = sample_ID + ' ' + father_ID + ' ' + mother_ID + '\n'
            o.write(out_line)
    f.close()
    o.close()

if __name__ == '__main__':
    main()
