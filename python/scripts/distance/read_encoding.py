def read_encoding_file(encoding_file):
    encodings = []
    f = open(encoding_file):
    for line in f:
        L = line.strip().split()
        encodings.append(L)
    f.close()
    return encodings
