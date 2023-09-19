from collections import deque
from collections import defaultdict
from collections import namedtuple
import networkx as nx
from matplotlib import pyplot as plt

ped_header = ['family_id',
              'individual_id',
              'paternal_id',
              'maternal_id',
              'sex',
              'phenotype']

Sample = namedtuple('Sample', ped_header)


def build_family_tree(ped_file):
    # ped_file_format: famID self father mother sex phenotype
    # read ped_file and build family tree
    G = {}
    with open(ped_file, 'r') as file:
        for line in file:
            values = line.strip().split()
            try:
                sample = Sample(*values)
            except TypeError:
                print('ERROR: PED file is not formatted correctly' + line)
            G[sample.individual_id] = [sample.paternal_id, sample.maternal_id]

    return G

def get_node_height(node, family_tree, family_members, root_p, root_m):
    # return height of node in G
    # height = number of edges from node to root
    # if node is root, height = 0
    height = 0
    queue = deque([node])
    visited = set()
    while queue:
        node = queue.popleft()
        if node not in visited:
            visited.add(node)
            neighbors = family_tree.get(node, [])
            for neighbor in neighbors:
                try:
                    neighbor_parents = family_tree.get(neighbor, [])
                    neighbor_father = neighbor_parents[0]
                except IndexError:
                    # neighbor is not in family tree
                    neighbor_father = None
                    pass
                # don't bother searching if neighbor is not in family (parent who married in)
                if neighbor not in visited and neighbor_father in family_members or neighbor == root_p or neighbor == root_m:
                    height += 1
                    queue.append(neighbor)
                    if neighbor == root_p or neighbor == root_m:
                        return height
    return height

def detect_unrelated(person, family_tree, root_p, root_m):
    # if there is a node that married in, their parents are not in the family
    # see if two neighbors of person are not in family_tree
    # if so, person is unrelated
    neighbors = family_tree.get(person, [])
    unrelated = False

    for neighbor in neighbors:
        if neighbor not in family_tree:
            if neighbor != root_p and neighbor != root_m:
                if person != root_p and person != root_m:
                    unrelated = True # person is unrelated
    return unrelated
