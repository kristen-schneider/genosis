import pysam

def read_data_file(data_file, delim=' '):
    data_dict = dict()
    f = open(data_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split(delim)
        seg = L[0]
        SNPs = int(L[1])
        data_dict[seg] = SNPs
    f.close()
    return data_dict

def vcf_samples(vcf_file):
    vcf_samples_list = []
    vcf_f = pysam.VariantFile(vcf_file)
    vcf_samples_list = list((vcf_f.header.samples))
    return vcf_samples_list
