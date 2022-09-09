import pysam

def read_endpoint_file(endpoint_file):
    '''
    Returns a dictionary whose key is a segment and
    whose value is a tuple with the start and endpoints
    basepbairs for that segment.
    '''

    seg_endpoint_dict = dict()

    f = open(endpoint_file, 'r')
    header = f.readline()
    for line in f:
        L = line.strip().split()
        seg = int(L[0])
        start_bp = int(L[1])
        end_pb = int(L[2])
        seg_endpoint_dict[seg] = (start_bp, end_pb)
    f.close()
    return seg_endpoint_dict

def vcf_samples(vcf_file):
    vcf_samples_list = []
    vcf_f = pysam.VariantFile(vcf_file)
    vcf_samples_list = list((vcf_f.header.samples))
    return vcf_samples_list


def sample_counts_from_vcf_block(vcf_file, chrm, start_bp, end_bp):
    '''
    Uses pysam to open vcf file with block compression.
    '''

    chrm = 'chr'+str(chrm)
    # gets sample list 
    vcf_samples_list = vcf_samples(vcf_file)
    
    sample_snp_dict = dict()

    # open vcf file and read by block record
    vcf_f = pysam.VariantFile(vcf_file)
    for snp in vcf_f.fetch(chrm, start_bp, end_bp):
        # get genotype for each sample
        for sample in vcf_samples_list:
            sample_GT = snp.samples[sample]['GT']
            if 1 in sample_GT:
                try:
                    sample_snp_dict[sample] += 1
                except KeyError:
                    sample_snp_dict[sample] = 1

    return sample_snp_dict
















