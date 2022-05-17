import matplotlib.pyplot as plt
import os

def plot_numElements(segcounts_dir, pltName):
    xy_dict = get_num_elements(segcounts_dir)

    plt.figure(figsize=(15,10))
    plt.plot(xy_dict.keys(), xy_dict.values(), 'o')
    plt.title('Number of elements by segment length')
    plt.xlabel('Segment Length')
    plt.ylabel('Total number of unique elements')
    plt.savefig(pltName)

def get_num_elements(segcounts_dir):
    xy_dict = dict()

    for f in os.listdir(segcounts_dir):
        f_open = open(segcounts_dir+f, 'r')
        num_elements = 0
        segment_length = 0
        for line in f_open:
            num_elements += 1
            A = line.strip().split()
            if segment_length == 0:
                segment = A[0]
                segment_length = len(segment)
        xy_dict[segment_length] = num_elements
    return xy_dict


