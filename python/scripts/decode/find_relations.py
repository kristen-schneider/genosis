from collections import deque
from collections import defaultdict
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

def search_family_graph(family_graph, family_members, i1, i2, family_roots):
    roots_i1 = family_roots[family_members[i1]]
    roots_i2 = family_roots[family_members[i2]]
    # if roots are different, return -1 (not in same family)
    if roots_i1 != roots_i2:
        return -1

    # search to find distance between individual_1 and individual_2
    # graph is of type FamilyNode

    # RETURN EARLY
    # if i1 or i2 not in any family, return None
    if i1 not in family_members.keys() or i2 not in family_members.keys():
        return None
    # if i1 or i2 not in the same family, return -1
    if family_members[i1] != family_members[i2]:
        return -1
    # if i1 == i2, return 0
    if i1 == i2:
        return 0

    # IN-LAWS
    # if current node's parent married in, keep their parents (otherwise removed downstream)
    for p in family_graph[i1].parents:
        # if parent is i2 return 1
        if p == i2:
            return 1
        # if parent is child of i2, return 2
        for c in family_graph[i2].children:
            if c == p:
                return 2
    # if current node married in, keep their children (otherwise removed downstream)
    for c in family_graph[i1].children:
        # if child is i2, return 1
        if c == i2:
            return 1
        # if child is parent of i2, return 2
        for p in family_graph[i2].parents:
            if p == c:
                return 2
    # if i1 married in (parent's parents are not in family), return -2
    if family_graph[i1].parents == []:
        return -2
    # if i2 married in (parent's parents are not in family), return -2
    if family_graph[i2].parents == []:
        return -2

    # BLOOD RELATIVES
    # find the shortest path between i1 and i2
    queue = deque()
    visited = {i: False for i in family_members.keys()}
    distances = {i: 0 for i in family_members.keys()}
    predecessors = {i: 0 for i in family_members.keys()}

    for i in family_members.keys():
        distances[i] = 1000000
        predecessors[i] = -1

    visited[i1] = True
    distances[i1] = 0
    queue.append(i1)

    # while the path has not been found
    while queue:
        current_node = queue.popleft()

        neighbors = get_neighbors(family_graph, family_graph[current_node], family_roots[family_members[current_node]], visited)
        for neighbor in neighbors:
            if not visited[neighbor]:
                visited[neighbor] = True
                distances[neighbor] = distances[current_node] + 1
                predecessors[neighbor] = current_node
                queue.append(neighbor)

                if neighbor == i2:
                    return distances[neighbor]

    # no path found (i1 and i2 are not related)
    return -1

def get_neighbors(family_graph, family_member, family_roots, visited):
    # return unvisited neighbors of family_member
    # do not include parents who married into family, or who do not have parents listed
    # do not include children who do not have parents listed
    parents_ = family_member.parents
    parents = [p for p in parents_]
    children_ = family_member.children
    children = [c for c in children_]

    if family_member.name in family_roots:
        parents = []
    else:
        for p in parents_:
            # remove parents that have already been visited
            if visited[p]:
                parents.remove(p)
                continue
            elif p in family_roots:
                continue
            else:
                if family_graph[p].parents == []:
                    parents.remove(p)
                    continue
    all_neighbors = parents + children
    return all_neighbors
