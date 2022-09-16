import argparse

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--plink')
    parser.add_argument('--samples_dir')
    parser.add_argument('--numQ')
    return parser.parse_args()

def main():
    args = get_args()

if __name__ == '__main__':
    main()
