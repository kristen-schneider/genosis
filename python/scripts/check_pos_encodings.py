import os
import argparse


MINIMUM=7   # minimum length for an encoding vector for the segment to be included

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pos_dir', type=str, help='directory of positional encoding vectors')
    parser.add_argument('--pos_ext', type=str, default='.pos', help='positional encoding vector extension')
    return parser.parse_args()

def main():
    args = parse_args()
    pos_dir = args.pos_dir
    pos_ext = args.pos_ext
    
    # look for files with positional encoding extension
    for pos_file in os.listdir(pos_dir):
        print(pos_file)
        if pos_file.endswith(pos_ext):
            pos_encoding = pos_dir + pos_file
            check_for_empty_entries(pos_encoding)

def check_for_empty_entries(pos_encoding):
    """
    check for segments with encoding vectors less than MINIMUM
    """
    #empty = False   # flag if positional encoding less than MINIMUM is found
    p = open(pos_encoding, 'r')
    for line in p:
        L = line.strip().split()
        # if len of positional encoding vector is less than MINIMUM:
        if len(L) < MINIMUM:
            print('!! Removing', pos_encoding)
            os.remove(pos_encoding) 
            #empty = True
            break       # exit loop as soon as bad encoding is found. don't keep searching file
    ## remove a segment file if there are bad encodings (less than MINIMUM)
    #if empty:
    #    print('!! Removing', pos_encoding)
    #    os.remove(pos_encoding) 
if __name__ == "__main__":
    main()
