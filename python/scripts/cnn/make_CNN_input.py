import sys

sample_IDs_file = sys.argv[1]
sample_encodings_file = sys.argv[2]
plink_IBD_file = sys.argv[3]
ID_CNN_file = sys.argv[4]
encoding_CNN_file = sys.argv[5]

def main():
    ID_list = sample_IDs_list(sample_IDs_file)
    ID_encoding_dictionary = ID_encoding_dict(ID_list, sample_encodings_file)
    write_CNN_input(ID_encoding_dictionary, plink_IBD_file, ID_CNN_file, encoding_CNN_file)

def sample_IDs_list(sample_IDs_file):
    ID_list = []
    
    ID_f = open(sample_IDs_file, 'r')
    for line in ID_f:
        ID_list.append(line.strip())

    ID_f.close()
    
    return ID_list

def ID_encoding_dict(ID_list, sample_encodings_file):
    
    ID_i = 0
    ID_encoding_dictionary = dict()
    encodings_f = open(sample_encodings_file, 'r')
   
    for line in encodings_f:
        encoding = line.strip()
        ID_encoding_dictionary[ID_list[ID_i]] = encoding
        ID_i += 1
    
    encodings_f.close()

    return ID_encoding_dictionary

def write_CNN_input(ID_encoding_dictionary, plink_IBD_file, ID_CNN_file, encoding_CNN_file):

    plink_f = open(plink_IBD_file, 'r')
    ID_header = 'sample1_ID\tsample2_ID\tdistance\n'
    ID_CNN_f = open(ID_CNN_file, 'w')
    ID_CNN_f.write(ID_header)
    encoding_header = 'sample1_encoding\tsample2_encoding\tdistance\n'
    encoding_CNN_f = open(encoding_CNN_file, 'w')
    encoding_CNN_f.write(encoding_header)

    header = None
    for line in plink_f:
        if (header == None):
            header = line
        else:
            A = line.strip().split()
            ID1 = A[1]
            ID2 = A[3]
            distance = float(A[11])

            ID1_encoding = ID_encoding_dictionary[ID1]
            ID2_encoding = ID_encoding_dictionary[ID2]
        
            ID_pairwise_entry = ID1 + '\t' + ID2 + '\t' + str(distance) + '\n'
            ID_CNN_f.write(ID_pairwise_entry)
            encoding_pairwise_entry = ID1_encoding + '\t' + ID2_encoding + '\t' + str(distance) + '\n'
            encoding_CNN_f.write(encoding_pairwise_entry)

    plink_f.close()
    ID_CNN_f.close()
    encoding_CNN_f.close()


if __name__ == '__main__':
    main()
