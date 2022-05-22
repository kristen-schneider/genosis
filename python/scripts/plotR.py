import matplotlib.pyplot as plt

def plot_frequencey(freq_counts_txt, pltName):
    f = read_freq_counts(freq_counts_txt)

    x = f.keys()
    y = f.values()

    plt.figure(figsize=(15, 10))
    # plt.yscale("log")
    # plt.xscale("log")
    # plt.bar(x,y, width=100)
    # plt.hist(y, bins=30)
    plt.bar(x, y)
    plt.title('Diversity of segments')
    plt.xlabel('Number of unique strings per segment')
    plt.ylabel('Frequency')
    plt.savefig(pltName)

def read_freq_counts(freq_counts_txt):
    header = None
    f = open(freq_counts_txt, 'r')
    if header == None:
        header = f.readline()

    freq_data_dict = {}
    for line in f:
        A = line.strip().split()
        numUniqueSegs = int(A[0])
        freq = int(A[1])
        freq_data_dict[numUniqueSegs] = freq

    return freq_data_dict