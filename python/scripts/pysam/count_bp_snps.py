import read_map

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--map')
    parser.add_argument('--cm-max')

    return parser.parse_args()

def main():
    args = get_args()
    


if __name__ == '__main__':
    main()
