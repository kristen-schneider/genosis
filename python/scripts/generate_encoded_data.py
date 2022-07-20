import random
import sys

numVariants = int(sys.argv[1])
numSamples = int(sys.argv[2])
outfile = sys.argv[3]

# does not represent true genomic data
def main():
    simulate_data(numVariants, numSamples, outfile)

def simulate_data(numVariants, numSamples, outfile):
    encodings = [0, 1, 2, 3]

    f = open(outfile, 'w')

    for s in range(numSamples):
        sample = random.choices(encodings, weights = [10, 5, 1, 1], k=numVariants)
        
        #for v in range(numVariants):
            #v = np.random.choice(np.arange(0, 3), p=[0.1, 0.05, 0.05, 0.2, 0.4, 0.2])
            #v = randint(0,3)
        #    sample += str(v)
        #print(sample)
        
        for v in sample:
            f.write(str(v))
        f.write('\n')
    
    f.close()


if __name__ == '__main__':
    main()
