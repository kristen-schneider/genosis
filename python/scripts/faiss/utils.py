def sum_distances(distances):
    penalty = 100
    sum = 0
    for d in distances:
        try:
            sum += d
        except TypeError:
            sum += penalty
    return sum
