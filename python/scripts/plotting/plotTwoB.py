import matplotlib.pyplot as plt

def plot_numElements(encoded_data, pltName):
    numVariants = get_num_variants(encoded_data)
    x = range(1, 3000)
    y = [(numVariants/i) for i in x]

    plt.figure(figsize=(15,10))
    plt.plot(x, y, 'o')
    plt.title('' + str(numVariants))
    plt.xlabel('Segment Length')
    plt.ylabel('Total number of elements')
    plt.savefig(pltName)

def get_num_variants(encoded_data):
    numVariants = 0
    f_open = open(encoded_data, 'r')
    oneLine = f_open.readline()
    for o in oneLine:
        if o == ' ':
            continue
        else: numVariants += 1
    return numVariants