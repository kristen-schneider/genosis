import math
import numpy as np
from scipy.spatial import distance
import os, sys, inspect

# euclidean distance
# hamming distance (hd and 1-hd)


def hamming_distance(v1, v2):
    '''
    counts the number of mismatches between 
    two input vectors
    '''
    hd = 0
    size_vector = len(v1)
    for v in range(size_vector):
        if v1[v] == v2[v]:
            continue
        else:
            hd += 1
    return hd



def euclidean_distance(v1, v2):
    '''
    computes the euclidean distance between 
    two input vectors
    '''
    size_vector = len(v1)
    ed = 0
    running_sum = 0

    for v in range(size_vector):
        diff = v1-v2
        diff_sqrd = pow(diff, 2)
        running_sum += diff_sqrd
    ed = math.sqrt(running_sum)
    return ed
