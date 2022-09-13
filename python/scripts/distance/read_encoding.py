def read_encoding_file(encoding_file):
    encodings = []
    f = open(encoding_file, 'r')
    for line in f:
        L = line.strip().split()
        sample_ID = L[0]
        encodings.append([float(i) for i in L[1:]])
    f.close()
    return encodings
