def write_seg_count(unique_seg_counts, seg_counts_txt):
    o = open(seg_counts_txt, 'w')
    for u_seg in unique_seg_counts.keys():
        s = u_seg + '\t' + str(unique_seg_counts[u_seg]) + '\n'
        o.write(s)

def get_data(seg_counts_txt):
    f = open(seg_counts_txt, 'r')
    data_dict = {}
    for line in f:
        A = line.strip().split()
        segment = A[0]
        count = int(A[1])
        data_dict[segment] = count

    return data_dict

def get_frequency(seg_counts_dict):
    freq_dict = {}
    for c in seg_counts_dict.values():
        try: freq_dict[c] += 1
        except KeyError: freq_dict[c] = 1

    return freq_dict

def write_freq(seg_counts_txt, freq_counts_txt):
    o = open(freq_counts_txt, 'w')
    seg_counts_dict = get_data(seg_counts_txt)
    freq_dict = get_frequency(seg_counts_dict)
    for count in freq_dict.keys():
        s = str(count) + '\t' + str(freq_dict[count]) + '\n'
        o.write(s)



