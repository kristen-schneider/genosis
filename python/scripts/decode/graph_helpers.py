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
            if sample.maternal_id in G:
                G[sample.maternal_id].children.append(sample.individual_id)
    return G

def search_family_graph(family_graph, family_members, i1, i2, root_p, root_m):
    # search to find distance between individual_1 and individual_2
    # graph is of type FamilyNode

    # if i1 or i2 not in family, return None
    if i1 not in family_members or i2 not in family_members:
        return None
    # if i1 == i2, return 0
    if i1 == i2:
        return 0

    # keep parents/children who are unrelated
    for p in family_graph[i1].parents:
        if p == i2:
            return 1
    for c in family_graph[i1].children:
        if c == i2:
            return 1

    # find shortest path between i1 and i2
    # if no path exists, return None
    queue = deque()
    visited = {i: False for i in family_members}
    distances = {i: 0 for i in family_members}
    predecessors = {i: 0 for i in family_members}

    for i in family_members:
        distances[i] = 1000000
        predecessors[i] = -1

    visited[i1] = True
    distances[i1] = 0
    queue.append(i1)

    while queue:
        current_node = queue.popleft()
        print(current_node)
        if family_graph[current_node].parents[0] not in family_members:
            if current_node == root_m or current_node == root_p:
                pass
            else:
                break
        neighbors = get_neighbors(family_graph, family_graph[current_node], root_p, root_m, visited)
        for neighbor in neighbors:
            if not visited[neighbor]:
                visited[neighbor] = True
                distances[neighbor] = distances[current_node] + 1
                predecessors[neighbor] = current_node
                queue.append(neighbor)

                if neighbor == i2:
                    return distances[neighbor]





def get_neighbors(family_graph, family_member, root_p, root_m, visited):
    # return neighbors of family_member

    parents_ = family_member.parents
    parents = [p for p in parents_]
    children_ = family_member.children
    children = [c for c in children_]

    # if person is a root node, remove their parents
    if family_member.name == root_p or family_member.name == root_m:
        parents = []
    else:
        for p in parents_:
            # remove parents that have already been visited
            if visited[p]:
                parents.remove(p)
                continue
            # don't include people who married into family
            elif p == root_p or p == root_m:
                continue
            else:
                try:
                    grandparents = family_graph[p].parents
                    if grandparents[0] not in family_graph or grandparents[1] not in family_graph:
                        parents.remove(p)
                except KeyError:
                    pass

    for c in children_:
        # remove children that have already been visited
        if visited[c]:
            children.remove(c)

    all_neighbors = parents + children
    return all_neighbors




def get_node_degree(node, G):
    # return degree of node in G
    # degree = number of edges connected to node
    degree = 0
    neighbors = G.get(node, [])
    degree += len(neighbors)
    return degree

def plot_graph(family_tree):
    # use networkx to plot family tree
    G = nx.DiGraph()
    for node in family_tree:
        G.add_node(node)
        for neighbor in family_tree[node]:
            G.add_edge(node, neighbor)

    # no overlap between nodes
    plt.figure(figsize=(20, 20))
    pos = nx.spring_layout(G, k=0.5, iterations=20)
    nx.draw(G, pos, with_labels=True)
    plt.savefig()
