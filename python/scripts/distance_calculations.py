import math
import numpy as np
from scipy.spatial import distance
import sys


encoding_file = sys.argv[1]
index1 = sys.argv[2]
index2 = sys.argv[3]

def main():
    
    print(brute_force(index1, index2))
    print(numpy_euclidean(index1, index2))

def brute_force(i1, i2):
    f = open(faiss_file, 'r')
    l_count = 0
    vector1, vector2 = '', ''

    for line in f:
        if l_count == i1:
            vector1 = line.strip()
        if l_count == i2:
            vector2 = line.strip()
        l_count += 1

    return euclidean_distance(vector1, vector2)

def numpy_euclidean(i1, i2):
    f = open(faiss_file, 'r')
    l_count = 0
    vector1, vector2 = '', ''
    vector1_list, vector2_list = [], []

    for line in f:
        if l_count == i1:
            vector1 = line.strip()
            for v in vector1: vector1_list.append(float(v))

        if l_count == i2:
            vector2 = line.strip()
            for v in vector2: vector2_list.append(float(v))

        l_count += 1

    # distance.euclidean(vector1_list, vector2_list)
    vector1_list = np.array(vector1_list)
    vector2_list = np.array(vector2_list)
    return np.sqrt(np.sum((vector1_list-vector2_list)**2))

def euclidean_distance(vector1, vector2):
    sum = 0
    for v in range(len(vector1)):
        v1 = int(vector1[v])
        v2 = int(vector2[v])

        diff = v1-v2
        diff_sqrd = pow(diff, 2)
        sum += diff_sqrd
    return math.sqrt(sum)

if __name__ == '__main__':
    main()

