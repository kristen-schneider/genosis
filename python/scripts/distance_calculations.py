import math
import numpy as np
from scipy.spatial import distance
import sys



def brute_force(i1, i2):
    f = open(encoding_file, 'r')
    l_count = 0
    vector1, vector2 = '', ''

    for line in f:
        if l_count == i1:
            vector1 = line.strip()
        if l_count == i2:
            vector2 = line.strip()
        l_count += 1
    
    f.close()
    return euclidean_distance(vector1, vector2)

def numpy_euclidean(i1, i2):
    f = open(encoding_file, 'r')
    l_count = 0
    vector1, vector2 = '', ''
    vector1_list, vector2_list = [], []

    for line in f:
        if l_count == i1:
            vector1_list = line.strip().split()
            if(len(vector1_list) == 1):
                vector1_list = []
                for i in range(len(line.strip())):
                    vector1_list.append(float(line[i]))
            

        if l_count == i2:
            vector2_list = line.strip().split()
            if(len(vector2_list) == 1):
                vector2_list = []
                for i in range(len(line.strip())):
                    vector2_list.append(float(line[i]))

        l_count += 1

    f.close()
    
    # distance.euclidean(vector1_list, vector2_list)
    vector1_list = np.array(vector1_list)
    vector2_list = np.array(vector2_list)
    return np.sqrt(np.sum((vector1_list-vector2_list)**2))

def euclidean_distance(vector1, vector2):
    sum = 0
    vector1_list = vector1.strip().split()
    vector2_list = vector2.strip().split()

    if(len(vector1_list) == 1):
        vector1_list = []
        vector2_list = []
        for i in range(len(vector1)):
            vector1_list.append(vector1[i])
            vector2_list.append(vector2[i])

    for v in range(len(vector1_list)):
        v1 = float(vector1_list[v])
        v2 = float(vector2_list[v])
        diff = v1-v2
        diff_sqrd = pow(diff, 2)
        sum += diff_sqrd
    return math.sqrt(sum)
