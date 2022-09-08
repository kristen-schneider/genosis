def searc_map_file(map_file=None, start_bp=0, end_bp=0):
    
    f = open(map_file, 'r')

    for line in f:
        L = line.strip().split()
        chrm = L[0]
        cm = float(L[2])
        bp = int(l[3])

    f.close()
