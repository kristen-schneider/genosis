import sys
from collections import defaultdict

#sys.path.insert(1, '~genotype-encoding/python/utils/')
#import distance_calculations
#from python.utils import distance_calculations

def get_sample_ID_list(sample_ID_file):
    all_sample_IDs = []
    f = open(sample_ID_file, 'r')
    for line in f:
        ID = line.strip()
        all_sample_IDs.append(ID)
    f.close()

    return all_sample_IDs

def get_encoding_list(sample_encodings_file):
    all_sample_encodings = []
    f = open(sample_encodings_file, 'r')
    for line in f:
        encoding = line.strip()
        encoding_list_ints = [int(i) for i in encoding]
        all_sample_encodings.append(encoding_list_ints)
    f.close()

    return all_sample_encodings

def get_ID_encoding_dict(sample_IDs, sample_encodings_file):
    ID_encoding_dict = dict.fromkeys(sample_IDs)

    ID_i = 0
    f = open(sample_encodings_file, 'r')
    for line in f:
        encoding_i = line.strip()
        ID_encoding_dict[sample_IDs[ID_i]] = encoding_i
        ID_i += 1
    f.close()

    if len(sample_IDs) != len(ID_encoding_dict.keys()):
        print('ERROR: not the same number of sample IDs and encodings.')
    return ID_encoding_dict

def get_query_euclidean_dict(query_ID, ID_encoding_dict):
    query_euclidean_dict = {ID: -1 for ID in ID_encoding_dict.keys()}
    query_encoding = ID_encoding_dict[query_ID]

    for ID in ID_encoding_dict.keys():
        sample_encoding = ID_encoding_dict[ID]
        euclidean_distance = distance_calculations.euclidean_distance(query_encoding, sample_encoding)
        query_euclidean_dict[ID] = euclidean_distance
    return query_euclidean_dict

def get_plink_dict(plink_file):
    plink_dict = defaultdict(dict)

    f = open(plink_file, 'r')
    header = f.readline()
    for line in f:
        l = line.strip().split()
        ID1 = l[1]
        ID2 = l[3]
        plink_distance = float(l[7])
        plink_dict[ID1][ID2] = plink_distance
        plink_dict[ID2][ID1] = plink_distance
    f.close()
    return plink_dict

def get_faiss_distances_dict(sample_ID_list, faiss_file):
    '''

    :param sample_ID_list:
    :param faiss_file:
    :return: return a dictionary of dictionaries
      whose first key is sample 1,
      second key is sample2 and
      whose value is the distance between the two samples
    '''
    faiss_dict = defaultdict(dict)
    num_samples = len(sample_ID_list)
    pairwise_count = 0

    f = open(faiss_file, 'r')
    for line in f:
        if "QUERY" in line:
            pairwise_count = 0
            sampleID1 = sample_ID_list[int(line.split(':')[1].strip())]
        else:
            if pairwise_count == num_samples:
                continue
            else:
                l = line.strip().split()
                try:
                    index = int(l[0])
                    sampleID2 = sample_ID_list[index]
                    distance = float(l[1])
                    faiss_dict[sampleID1][sampleID2] = distance
                    pairwise_count += 1

                except ValueError:
                    continue
    f.close()
    return faiss_dict

def get_faiss_rankings(sample_ID_list, faiss_file):
    '''

    :param sample_ID_list:
    :param faiss_file:
    :return: dictionary whose key is a sample ID
    and whose value is the ordered list of rankings
    '''
    num_samples = len(sample_ID_list)
    faiss_rankings = {sampleID: [] for sampleID in sample_ID_list}
    f = open(faiss_file, 'r')

    value_ranked_indexes = []
    for line in f:
        if "QUERY" in line:
            sampleID = sample_ID_list[int(line.split(':')[1].strip())]
        else:
            if len(value_ranked_indexes) == num_samples:
                faiss_rankings[sampleID] = value_ranked_indexes
                value_ranked_indexes = []
            else:
                l = line.strip().split()
                try:
                    index = int(l[0])
                    value_ranked_indexes.append(index)
                except ValueError:
                    continue
    f.close()
    return faiss_rankings


