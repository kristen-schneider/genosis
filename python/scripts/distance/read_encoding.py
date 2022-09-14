def read_encoding_file(encoding_file):
    encodings = dict()
    f = open(encoding_file, 'r')
    for line in f:
        print(line)
        L = line.strip().split()
        sample_ID = L[0]
        encodings[sample_ID] = [float(i) for i in L[1:]]
    f.close()
    return encodings
