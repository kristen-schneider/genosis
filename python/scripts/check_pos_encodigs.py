import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pos_dir', type=str)
    parser.add_argument('--pos_ext', type=str)
    return parser.parse_args()

def main():
    args = parse_args()
    pos_dir = args.pos_dir
    pos_ext = args.pos_ext

    for pos_file in os.listdir(pos_dir):
        if pos_file.endswith(pos_ext):
            #print(pos_file)
            pos_encoding = pos_dir + pos_file
            check_for_empty_entries(pos_encoding)

def check_for_empty_entries(pos_encoding):
    empty = False
    p = open(pos_encoding, 'r')
    for line in p:
        L = line.strip().split()
        if len(L) < 5:
            empty = True
            print(line)
            break
    if empty:
        print('!! Removing ', pos_encoding, ' for empty positional encodings')
        os.remove(pos_encoding) 

if __name__ == "__main__":
    main()
