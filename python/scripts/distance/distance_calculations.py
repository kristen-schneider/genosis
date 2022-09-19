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
        diff = v1[v] - v2[v]
        diff_sqrd = pow(diff, 2)
        running_sum += diff_sqrd
    ed = math.sqrt(running_sum)
    return ed

def kristen(v1, v2, gaps_allowed):
    kd = 0
    running_score = 0
    num_consecutive_matches = 0
    gaps = 0
    all_scores = []

    size_vector = len(v1)
    for v in range(size_vector):
        if int(v1[v]) == 1 and int(v1[v]) == int(v2[v]):
            num_consecutive_matches += 1
            running_score += num_consecutive_matches
        else:
            gaps += 1
            if gaps > gaps_allowed:
                if running_score != 0:
                    all_scores.append(running_score)
                running_score = 0
                num_consecutive_matches = 0
                gaps = 0
            else:
                continue
    if running_score != 0:
        all_scores.append(running_score)
    kd = sum(all_scores)
    return kd
