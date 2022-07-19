vcf_file = '/home/sdp/precision-medicine/data/vcf/ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf'
encoded_file = '/home/sdp/precision-medicine/data/encoded/new.encoded.txt'

def main():
    samples = sample_name_index(vcf_file)
    dict_samples = sample_encoding_dict(samples, encoded_file)


def sample_encoding_dict(samples, encoded_file):
    f = open(encoded_file, 'r')
    dict_samples = dict()

    line_i = 0
    for line in encoded_file:
        sampleID = samples[line_i]
        dict_samples[sampleID] = line
    
    return dict_samples


def sample_name_index(vcf_file):
    f = open(vcf_file, 'r')
    samples = []

    for line in f:
        if line.startswith('##'):
            continue
        else:
            if line.startswith('#CHROM'):
                A = line.strip().split()
                samples = [9:]
                return samples

if __name__ == '__main__':
    main()


