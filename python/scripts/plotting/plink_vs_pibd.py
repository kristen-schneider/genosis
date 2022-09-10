import argparse
import matplotlib.pyplot as plt

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pibd_plink_file')
    parser.add_argument('--out_png')
    parser.add_argument('--chrm')
    parser.add_argument('--method')
    return parser.parse_args()

def main():
    args = get_args()
    [plink_x, pibd_y] = read_file(args.pibd_plink_file)
    plot_pibd_plink(plink_x, pibd_y, args.chrm, args.out_png, args.method)


def read_file(pibd_plink_file):
    plink_x = []
    pibd_y = []
    f = open(pibd_plink_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split()
        ilash = float(L[0])
        plink = float(L[1])
        if plink == -1:
            continue
        plink_x.append(plink)
        pibd_y.append(ilash)

    f.close()
    return [plink_x, pibd_y]

def plot_pibd_plink(plink_x, pibd_y, chrm, out_png, method):
    if method == ('agg'):
        png_name = out_png + 'chr' + str(chrm) \
                   + '.plink-pibd.agg.png'
        title = 'Plink vs. Phased IBD for full VCF\n' \
                '(Chromosome ' + str(chrm) + ')'
    elif method == ('ind'):
        png_name = out_png + 'chr' + str(chrm) \
                   + '.plink-pibd.ind.png'
        title = 'Plink vs. Phased IBD for IBD segments reported by Phased IBD\n' \
                '(Chromosome ' + str(chrm) + ')'

    plt.figure()#figsize=(15, 12))
    ax1 = plt.subplot(111)
    ax1.set_title(title)
    ax1.scatter(plink_x, pibd_y, color='steelblue')
    ax1.set_xlabel('Plink')
    # ax1.set_xlim([0, 1])
    ax1.set_ylabel('Phased IBD\n(cM distance)')
    # ax1.set_ylim([0, max(ilash_y) + 5])
    # ax1.set_yscale('log')
    ax1.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    plt.savefig(png_name)

if __name__ == '__main__':
    main()
