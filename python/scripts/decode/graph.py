from collections import namedtuple

ped_header = ['family_id',
              'individual_id',
              'paternal_id',
              'maternal_id',
              'sex',
              'phenotype']

Sample = namedtuple('Sample', ped_header)

class FamilyNode:
    # has parents and children
    def __init__(self, name):
        self.name = name
        self.parents = []
        self.children = []

def build_family_graph(ped_file):
    # ped_file_format: famID self father mother sex phenotype
    # read ped_file and build family graph (bidirectional)
    # create a graph of type FamilyNode
    G = {}
    with open(ped_file, 'r') as file:
        for line in file:
            values = line.strip().split()
            try:
                sample = Sample(*values)
            except TypeError:
                print('ERROR: PED file is not formatted correctly' + line)

            # create node from this sample
            sample_node = FamilyNode(sample.individual_id)
            # add parents to node
            sample_node.parents.append(sample.paternal_id)
            sample_node.parents.append(sample.maternal_id)
            # add node to graph
            G[sample.individual_id] = sample_node

            # add children to parents
            if sample.paternal_id in G:
                G[sample.paternal_id].children.append(sample.individual_id)
            else:
                sample_node = FamilyNode(sample.paternal_id)
                sample_node.children.append(sample.individual_id)
                G[sample.paternal_id] = sample_node
            if sample.maternal_id in G:
                G[sample.maternal_id].children.append(sample.individual_id)
            else:
                sample_node = FamilyNode(sample.maternal_id)
                sample_node.children.append(sample.individual_id)
                G[sample.maternal_id] = sample_node
    return G

def get_distance_to_root(family_graph, i1, root, distance):
    # if no parents, return -1
    if len(family_graph[i1].parents) == 0:
        return -1
    if i1 == root:
        return distance
    else:
        distance += 1
        curr_parents = family_graph[i1].parents

        for p in curr_parents:
            if family_graph[p].parents != []:
                germline_parent = p
                if germline_parent == root:
                    return distance
            # else:
            #     continue
        distance = get_distance_to_root(family_graph, germline_parent, root, distance)
        # for parent in family_graph[i1].parents:
        #     get_distance_to_root(family_graph, parent, root, distance)
    return distance

def get_family_members(ped_file):
    family_members = dict()
    with open(ped_file, 'r') as file:
        for line in file:
            values = line.strip().split()
            family_ID = values[0]
            samples = values[1:4]
            for sample in samples:
                family_members[sample] = family_ID
    return family_members

def get_family_roots(roots_file):
    family_roots = dict()
    with open(roots_file, 'r') as file:
        for line in file:
            values = line.strip().split()
            family_ID = values[0]
            root_node = values[1]
            family_roots[family_ID] = root_node
    return family_roots

