import argparse
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ilash_plink_file')
    parser.add_argument('--out_png')
    parser.add_argument('--chrm')
    parser.add_argument('--method')
    return parser.parse_args()

def main():
    args = get_args()
    [plink_x, ilash_y] = read_file(args.ilash_plink_file)
    plot_ilash_plink(plink_x, ilash_y, args.chrm, args.out_png, args.method)


def read_file(ilash_plink_file):
    plink_x = []
    ilash_y = []
    f = open(ilash_plink_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split()
        ilash = float(L[0])
        plink = float(L[1])
        if plink == -1:
            continue
        plink_x.append(plink)
        ilash_y.append(ilash)

    f.close()
    return [plink_x, ilash_y]

def plot_ilash_plink(plink_x, ilash_y, chrm, out_png, method):
    if method == ('agg'):
        png_name = out_png + 'chr' + str(chrm) \
                   + '.plink-ilash.agg.png'
        title = 'Plink vs. iLASH for full VCF\n' \
                '(Chromosome ' + str(chrm) + ')'
    elif method == ('ind'):
        png_name = out_png + 'chr' + str(chrm) \
                   + '.plink-ilash.ind.png'
        title = 'Plink vs. iLASH for IBD segments reported by iLASH\n' \
                '(Chromosome ' + str(chrm) + ')'

    plt.figure()#figsize=(15, 12))
    ax1 = plt.subplot(111)
    ax1.set_title(title)
    ax1.scatter(plink_x, ilash_y, color='steelblue')
    ax1.set_xlabel('Plink')
    # ax1.set_xlim([0, 1])
    ax1.set_ylabel('iLASH\n(cM distance)')
    # ax1.set_ylim([0, max(ilash_y) + 5])
    # ax1.set_yscale('log')
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    plt.savefig(png_name)

if __name__ == '__main__':
    main()
