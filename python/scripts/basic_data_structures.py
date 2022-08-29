import sys
import re
from collections import defaultdict

# sys.path.insert(1, '~genotype-encoding/python/utils/')
# import distance_calculations
from python.utils import distance_calculations

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
        try:
            plink_distance = float(l[9])
        except ValueError:
            plink_distance = -1
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

def get_faiss_rankings_dict(sample_ID_list, faiss_file):
    '''

    :param sample_ID_list:
    :param faiss_file:
    :return: return a dictionary of dictionaries
      whose first key is sample 1,
      second key is sample2 and
      whose value is the rank for sample 2
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
                    rank = pairwise_count
                    sampleID2 = sample_ID_list[index]
                    faiss_dict[sampleID1][sampleID2] = rank
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


def get_plink_distances(database_IDs, queryIDs, plink_f):
    header_len = 1
    P = {}
    with open(plink_f) as f:
        for l in f:
            if "time" in l:
                break;
            if header_len > 0:
                header_len = header_len - 1
                continue
            A = l.rstrip().split()
            a = A[1]
            b = A[3]
            dist = float(A[11])

            if a in queryIDs:
                if a not in P:
                    P[a] = []
                if b in database_IDs:
                    P[a].append((b, dist))

            if b in queryIDs:
                if b not in P:
                    P[b] = []
                if a in database_IDs:
                    P[b].append((a, dist))

        f.close()
        # sort by distance
        for p in P:
            P[p].sort(key=lambda tup: tup[1], reverse=True)

        return P

def get_faiss_distances(database_IDs, queryIDs, faiss_file):
    F = {}
    query_id = None
    f = open(faiss_file, 'r')
    for l in f:
        if len(l) == 1:
            continue
        elif re.search(r'^TIME', l):
            continue
        elif re.search(r'^QUERY', l):
            query_id = l.rstrip().split()[1]
            # query_id = int(l.rstrip().split()[1])
        elif query_id is not None and query_id in queryIDs:
            A = l.rstrip().split()
            idx_ID = A[0]
            idx = int(A[1])
            dist = float(A[2])
            if query_id not in F:
                F[query_id] = []
            # if len(F[query_id]) > 50:
            #     query_id = None
            #     continue
            # if idx_ID in database_IDs:
            F[query_id].append((idx_ID, dist))
    f.close()

    # sort by distance
    for f in F:
        F[f].sort(key=lambda tup: tup[1])
    return F

def get_db_q_IDs(file):
    IDs = []
    # getting IDs of database samples
    with open(file) as f:
        for l in f:
            IDs.append(l.rstrip())
    return IDs

