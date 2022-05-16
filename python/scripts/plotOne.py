import matplotlib.pyplot as plt

def plot_frequencey(freq_counts_txt, pltName):
    f = read_freq_counts(freq_counts_txt)

    x = f.keys()
    y = f.values()


    plt.figure(figsize=(15,10))
    plt.yscale("log")
    plt.xscale("log")
    # plt.bar(x,y, width=100)
    # plt.hist(y, bins=30)
    plt.plot(x, y, 'o')
    plt.title('Diversity of segments ')
    plt.xlabel('Number of elements per bin\n(bin = segment)')
    plt.ylabel('Frequency')
    plt.savefig(pltName)

def read_freq_counts(freq_counts_txt):
    f = open(freq_counts_txt, 'r')
    freq_data_dict = {}
    for line in f:
        A = line.strip().split()
        count = int(A[0])
        freq = int(A[1])
        freq_data_dict[count] = freq

    return freq_data_dict

