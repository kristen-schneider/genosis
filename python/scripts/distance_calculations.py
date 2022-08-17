import math

def euclidean_distance(vector1, vector2):
    sum = 0
    for v in range(len(vector1)):
        v1 = int(vector1[v])
        v2 = int(vector2[v])

        diff = v1-v2
        diff_sqrd = pow(diff, 2)
        sum += diff_sqrd
    return math.sqrt(sum)

