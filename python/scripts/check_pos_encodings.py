import os
import argparse

MINIMUM=7   # minimum length for an encoding vector for the segment to be included

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pos_enc', type=str, help='positional encoding file')
    return parser.parse_args()

def main():
    args = parse_args()
    pos_encoding_file = args.pos_enc
    check_for_empty_entries(pos_encoding_file)

def check_for_empty_entries(pos_encoding_file):
    """
    check for segments with encoding vectors less than MINIMUM
    """
    #empty = False   # flag if positional encoding less than MINIMUM is found
    p = open(pos_encoding_file, 'r')
    print(pos_encoding_file)
    for line in p:
        L = line.strip().split()
        # if len of positional encoding vector is less than MINIMUM:
        if len(L) < MINIMUM:
            print('!! Removing', pos_encoding_file)
            os.remove(pos_encoding_file) 
            break       # exit loop as soon as bad encoding is found. don't keep searching file

if __name__ == "__main__":
    main()
