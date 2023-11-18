import sys
import glob
import concurrent.futures
import os
import argparse
import numpy as np

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--id_file', type=str, required=True)
    parser.add_argument('--data_dir', type=str, required=True)
    parser.add_argument('--num_threads', type=int, default=1)
    parser.add_argument('--top_N', type=int, default=10)
    parser.add_argument('--out_file', type=str, required=True)
    return parser.parse_args()


def get_name_id_maps(id_file):
    name_id_map = {}
    id_name_map = {}
    with open(id_file) as lines:
        i = 0
        for line in lines:
            sample_name = line.rstrip()[:-2] # strip off happlotype 
            if sample_name not in name_id_map:
                name_id_map[sample_name] = i
                id_name_map[i] = sample_name
                i+=1
    return name_id_map, id_name_map

def process_file(file_path, name_id_map):
    sys.stderr.write(f'Processing file {file_path}\n')
    pop_counts = []
    for i in range(len(name_id_map)):
        pop_counts.append( np.array([0] * len(name_id_map) ))

    with open(file_path, 'r') as lines:
        target = None
        for line in lines:
            A = line.rstrip().split()
            if len(A) == 0 : continue
            if A[0] == 'Query:':
                target = A[1][:-2] # strip off happlotype 
            else:
                hit = A[0][:-2] # strip off happlotype 
                sim_score = float(A[1])
                
                pop_counts[name_id_map[target]][name_id_map[hit]] =  \
                        pop_counts[name_id_map[target]][name_id_map[hit]] + 1
    return pop_counts

def process_files_in_directory(data_dir, name_id_map, num_threads):
    file_paths = [os.path.join(data_dir, filename)
                  for filename in os.listdir(data_dir)
                  if os.path.isfile(os.path.join(data_dir, filename))]

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) \
            as executor:
        future_to_file = \
            {executor.submit(process_file, file_path, name_id_map): file_path \
            for file_path in file_paths}

        for future in concurrent.futures.as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                data = future.result()
                results.append(data)  # Aggregate results
            except Exception as exc:
                sys.stderr.write(f'{file_path} generated an exception: {exc}\n')

    return results

def get_top(counts, id_name_map, N):
    indexed_list = sorted(enumerate(counts), key=lambda x: x[1], reverse=True)
    top_N_indexes = [index for index, value in indexed_list[:N]]

    return [(id_name_map[i], counts[i])  for i in top_N_indexes]


def process_sample(i, id_name_map, all_hits, top_N):
    target_sample = id_name_map[i]
    sys.stderr.write(f'Processing sample {target_sample}\n')
    r = np.zeros(len(id_name_map))

    for hit in all_hits:
        r += hit[i]

    top_hits = get_top( r, id_name_map, top_N )
    return([target_sample, top_hits])

def process_samples(all_hits, id_name_map, num_threads, top_N):
    ids = range(len(id_name_map))
    results = [] 

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_threads) \
            as executor:
        future_to_file = \
            {executor.submit(process_sample,
                             i,
                             id_name_map,
                             all_hits,
                             top_N): \
             i for i in ids}

        for future in concurrent.futures.as_completed(future_to_file):
            i = future_to_file[future]
            try:
                data = future.result()
                results.append(data)  # Aggregate results
            except Exception as exc:
                sys.stderr.write(f'{i} generated an exception: {exc}')

    return results

def main():
    args = get_args()
    name_id_map, id_name_map = get_name_id_maps(args.id_file)
    all_hits = process_files_in_directory(args.data_dir,
                                          name_id_map,
                                          args.num_threads)

    all_tops = process_samples(all_hits,
                               id_name_map,
                               args.num_threads,
                               args.top_N)

    with open(args.out_file, 'w') as f:
        for top in all_tops:
            f.write('\t'.join([top[0]] \
                    + [','.join(str(y) for y in x) for x in top[1]]) \
                    + '\n')
                             

if __name__ == '__main__':
    main()
