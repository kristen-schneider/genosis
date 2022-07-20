vcf_file = '/home/sdp/precision-medicine/data/vcf/ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf'
encoded_file = '/home/sdp/precision-medicine/data/encoded/new.encoded.txt'
IBD_file = '/home/sdp/precision-medicine/data/IBD/ALL.chr14.genome'
cnn_file = '/home/sdp/precision-medicine/data/CNN.input.small.txt'

def main():
    samples = sample_name_index(vcf_file)
    samples_encoding = sample_encoding_dict(samples, encoded_file)

    samples_IBD = samples_IBD_dict(IBD_file, samples_encoding, cnn_file)

def samples_IBD_dict(IBD_file, samples_encoding, cnn_file):
    f = open(IBD_file, 'r')
    o = open(cnn_file, 'w')
    o.write('sample1\tsample2\tdistance\n')

    for line in f:
        A = line.strip().split()
        sample1 = A[1]
        sample2 = A[3]
        distance = A[11]
        try:
            sample1_encoding = samples_encoding[sample1]
            sample2_encoding = samples_encoding[sample2]
            #new_string = sample1_encoding.strip() + '\t'\
            #             + sample2_encoding.strip() + '\t'\
            #             + distance.strip() + '\n'
            
            new_string = sample1 + '\t' + sample2 + '\t' + distance + '\n'
            o.write(new_string)
        except KeyError:
            continue



def sample_encoding_dict(samples, encoded_file):
    f = open(encoded_file, 'r')
    dict_samples = dict()

    line_i = 0
    for line in f:
        sampleID = samples[line_i]
        dict_samples[sampleID] = line
        line_i += 1

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
                samples = A[9:]
                return samples

if __name__ == '__main__':
    main()

