from collections import namedtuple
from collections import deque
import argparse

ped_header = ['family_id',
              'individual_id',
              'paternal_id',
              'maternal_id',
              'sex',
              'phenotype']

Sample = namedtuple('Sample', ped_header)

def parse_args():
    parser = argparse.ArgumentParser(description='Get relations between individuals')
    parser.add_argument('--ped', type=str, help='PED file',
                        default='/Users/kristen/PycharmProjects/data_analyses/input_data/AFR-chrm10-everyone.fam')
    return parser.parse_args()

def main():
    args = parse_args()
    ped_file = args.ped
    relations_file = ped_file.replace('.fam', '.relations')

    samples = read_file(ped_file)
    G = make_graph(samples)
    V = G.keys()
    for v_1 in V:
        for v_2 in V:
            if v_1 != v_2:
                T_1, d_1 = bfs_tree(v_1, G)
                T_2, d_2 = bfs_tree(v_2, G)

                common = []
                for t in T_1:
                    if t in T_2:
                        common.append(t)

                depths = []
                for v in common:
                    depths.append((d_1[v], d_2[v]))
                if len(depths) > 0:
                    rce_depths = min(depths, key=lambda x: sum(x))
                    # print(rce_depths, name_rce(rce_depths))
                    r_depths, relation = name_rce(rce_depths)
                    direction = (v_1, v_2) if rce_depths == r_depths else (v_2, v_1)
                    # write to file
                    with open(relations_file, 'a') as file:
                        file.write(direction[0] + ',' + direction[1] + ',' + relation + '\n')
                    print(direction[0], direction[1], relation)


def read_file(filename):
    samples = []
    with open(filename, 'r') as file:
        for line in file:
            values = line.strip().split()
            sample = Sample(*values)
            samples.append(sample)
    return samples

def make_graph(samples):
    G = {}
    for s in samples:
        if s.paternal_id != '0' and s.maternal_id != '0':
            G[s.individual_id] = [s.paternal_id,s.maternal_id]
    return G

def bfs_tree(v, G):
    visited = set()
    queue = deque([v])

    T = {}
    depth = {v:0}

    while queue:
        node = queue.popleft()
        if node not in visited:
            visited.add(node)
            neighbors = G.get(node, [])
            for neighbor in neighbors:
                if neighbor not in visited:
                    if node not in T:
                        T[node] = []
                    T[node].append(neighbor)
                    depth[neighbor] = depth[node] + 1
                    queue.append(neighbor)
    return T, depth

def name_rce(depths):
    q = depths[0]
    t = depths[1]

    if q == t: #same generation
        if q == 1:
            return (depths, 'Sibling')
        else:
            return (depths, str(q-1) + ' Cousin')
    if q > t:  #earlier generatin
        if t == 0:
            grands = 1 if q - t >= 2 else 0
            greats = max(0, q - 3)
            return (depths, str(greats) + 'x  Great ' + str(grands) + 'x Grand Parent')
        if t == 1:
            grands = 1 if q - t >= 2 else 0
            greats = max(0, q - 3)
            return (depths, str(greats) + 'x  Great ' + str(grands) + 'x Grand Aunt_Uncle')
        else:
            cousin = t - 1
            removed = q - t
            return (depths, str(cousin) + 'x  Cousin ' + str(removed) + 'x Removed')
    if q < t:  #later generatin
        return name_rce((depths[1], depths[0]))

if __name__ == '__main__':
    main()

