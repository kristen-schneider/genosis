import matplotlib.pyplot as plt
import sys

sys.path.insert(1, '/home/sdp/genotype-encoding/python/utils/')
import basic_datastructures

sample_ID_file = sys.argv[1]
sample_encodings_file = sys.argv[2]
plink_file = sys.argv[3]
plink_euclidean_png = sys.argv[4]

def main():
    genotype_plink()

def genotype_plink():
    sampleIDs = basic_datastructures.get_sample_ID_list(sample_ID_file)
    ID_encoding_dict = basic_datastructures.get_ID_encoding_dict(sampleIDs, sample_encodings_file)

    query_ID = sampleIDs[0]
    query_encoding = ID_encoding_dict[query_ID]
    query_euclidean_dict = basic_datastructures.get_query_euclidean_dict(query_ID, ID_encoding_dict)
    query_plink_dict = basic_datastructures.get_plink_dict(plink_file)[query_ID]

    euclidean_data = []
    plink_data = []
    for sample in sampleIDs:
        # plink doesnt report self comparison
        try:
            plink_data.append(query_plink_dict[sample])
            euclidean_data.append(query_euclidean_dict[sample])
        except KeyError:
            continue

    #plt.figure(figsize=(20,20))
    plt.scatter(plink_data, euclidean_data)
    plt.title('PLINK IBD vs EUCLIDEAN DISTANCE')
    plt.xlabel('PLINK SIMILARITY')
    plt.ylabel('EUCLIDEAN DISTANCE')
    plt.savefig(plink_euclidean_png)

if __name__ == '__main__':
    main()
