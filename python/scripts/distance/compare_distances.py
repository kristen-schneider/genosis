import argparse
import distance_calculations

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--encoding_file')
    return parser.parse_args()

def main():
    args = get_args()


if __name__ == '__main__':
    main()
