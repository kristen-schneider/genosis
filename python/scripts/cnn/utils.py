
def get_dimensions(sample_encodings_file):
    f = open(sample_encodings_file, 'r')
    num_samples = 0
    for line in f:
        encoding = line.strip()
        num_variants = len(encoding)
        num_samples += 1
    f.close()
    return [num_variants, num_samples]

def count_pairs(num_samples):
    pairs = 0
    for i in range(num_samples):
        pairs += i
    return pairs

# def CNN_input_dimensions(tf_record):
#     """
#     find dimensions of incoming data
#         number of variants = size of input vector
#         number of lines = number of distances computed
#     """
#     num_variants = -1
#     num_lines = -1
#
#     f = open(tf_record, 'r')
#     header = f.readline()
#     num_lines = 0
#     for line in f:
#         if num_lines == 0:
#             d = line.strip().split()
#             sample1_encoded = d[0]
#             sample2_encoded = d[1]
#             if len(sample1_encoded) != len(sample2_encoded):
#                 print('Error. Two encodings of differnt length.')
#                 return -1
#             else:
#                 num_variants = len(sample1_encoded)
#         num_lines += 1
#     f.close()
#     return num_variants, num_lines
