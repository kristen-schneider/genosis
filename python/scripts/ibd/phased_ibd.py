from IPython.display import display
import phasedibd as ibd
import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--map')
    return parser.parse_args()

def main():
    args = get_args()
    
    print('vcf: ', args.vcf)
    print('map: ', args.map)
    
    #haplotypes = ibd.VcfHaplotypeAlignment(args.vcf)
    haplotypes = ibd.VcfHaplotypeAlignment(args.vcf, args.map)
    tpbwt = ibd.TPBWTAnalysis()
    ibd_results = tpbwt.compute_ibd(haplotypes)
    display(ibd_results)

if __name__ == '__main__':
    main()
