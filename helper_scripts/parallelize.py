import glob
from os.path import basename

vcf_dir='/Users/krsc0813/precision-medicine/example/vcf_segments/'

vcf_segments = glob.glob(vcf_dir + "*.vcf.gz")
vcf_segments = list(map(basename, vcf_segments))
vcf_segments = [".".join(v.split('.')[:-2]) for v in vcf_segments]

print(list(vcf_segments))
