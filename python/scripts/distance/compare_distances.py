import argparse
import distance_calculations
import read_encoding_file

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoding_file')
    return parser.parse_args()

def main():
    args = get_args()
    encodings = read_encoding_file.read_encoding_file(args.encoding_file)

    Query = encodings[0]
    Database = encodings


if __name__ == '__main__':
    main()
