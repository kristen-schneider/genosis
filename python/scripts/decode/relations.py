import argparse
from collections import defaultdict
import graph as g

def parse_args():
    parser = argparse.ArgumentParser(description='Label pedigree relations between individuals')
    parser.add_argument('--ped', type=str, help='PED (fam) file')
    parser.add_argument('--root_file', type=str, help='root nodes for family members')
    parser.add_argument('--relations', type=str, help='output relations file')
    parser.add_argument('--gen', type=int, help='number of generations', default=5)
    return parser.parse_args()

def main():
    args = parse_args()
    ped_file = args.ped
    relations_file = args.relations
    root_file = args.root_file
    gen = args.gen

    relations_labels = defaultdict(dict)

    family_graph = g.build_family_graph(ped_file)
    family_members = g.get_family_members(ped_file)
    family_root_nodes = g.get_family_roots(root_file)

    for i1 in family_members:
        family_i1 = family_members[i1]
        root_i1 = family_root_nodes[family_i1]
        for i2 in family_members:
            family_i2 = family_members[i2]
            root_i2 = family_root_nodes[family_i2]

            # check if same family, if not, return -1
            if root_i1 != root_i2:
                relation = 'unrelated'
            else:
                root_node = root_i1

                height_i1 = g.get_distance_to_root(family_graph, i1, root_node, distance=0)
                height_i2 = g.get_distance_to_root(family_graph, i2, root_node, distance=0)
                relation = label_relations(family_graph, i1, i2, height_i1, height_i2, gen)

            print(i1, i2, height_i1, height_i2, relation)
            try:
                relations_labels[i1][i2] = relation
            except KeyError:
                relations_labels[i1] = {i2: relation}

    # write to file
    with open(relations_file, 'w') as file:
        file.write('---------------------\n')
        file.write('i2 is the _____ of i1\n')
        file.write('---------------------\n')
        for i1 in family_members:
            for i2 in family_members:
                file.write(family_members[i1] + ',' + i1 + ',' + i2 + ',' + relations_labels[i1][i2] + '\n')

    return 0

def label_relations(family_graph, i1, i2, height_i1, height_i2, gen):
    relationship = 'undetermined'

    diff = height_i1 - height_i2

    # if diff > 1 grand = 1, else grand = 0
    if abs(diff) > 1:
        grand = 1
    else:
        grand = 0
    great = abs(diff) - 1
    cousin = 0
    removed = abs(diff)

    if diff == 0:
        my_parental_gen = get_parents(family_graph, [i1], abs(diff))
        relative_parental_gen = get_parents(family_graph, [i2], abs(diff))
        # self
        if i1 == i2:
            relationship = 'self'
            return relationship
        # siblings
        elif check_shared_parents(my_parental_gen, relative_parental_gen):
            relationship = 'sibling'
            return relationship
        # x-cousin
        for g in range(height_i1):
            cousin += 1
            my_parental_gen = get_parents(family_graph, my_parental_gen, 1)
            relative_parental_gen = get_parents(family_graph, relative_parental_gen, 1)
            if check_shared_parents(my_parental_gen, relative_parental_gen):
                relationship = str(cousin) + '-cousin'
                return relationship
        else:
            relationship = 'undetermined'

    # i1 is an ancestor of i2
    elif diff <= 0:
        my_parental_gen = get_parents(family_graph, [i1], abs(diff)-1)
        relative_parental_gen = get_parents(family_graph, [i2], abs(diff))

        for i in range(abs(diff)):
            # direct x-great x-grand child
            if i1 in relative_parental_gen:
                relationship = (great-1) * 'great-' + (grand) * 'grand' + 'child'
                return relationship
            # direct niece-nephew
            # increase relatives parents by 1
            relative_parental_gen = get_parents(family_graph, relative_parental_gen, 1)
            if check_shared_parents(my_parental_gen, relative_parental_gen):
                relationship = (great) * 'great-' + 'niece-nephew'
            # some kind of cousin
            else:
                for g in range(height_i1):
                    cousin += 1
                    my_parental_gen = get_parents(family_graph, my_parental_gen, 1)
                    relative_parental_gen = get_parents(family_graph, relative_parental_gen, 1)
                    if check_shared_parents(my_parental_gen, relative_parental_gen):
                        relationship = str(cousin) + '-cousin-' + str(removed) + '-removed'
                        return relationship
            grand = 1
            great += 1

    # i1 is a descendant of i2
    elif diff >= 0:
        my_parental_gen = get_parents(family_graph, [i1], abs(diff))
        relative_parental_gen = get_parents(family_graph, [i2], abs(diff)-1)

        for i in range(abs(diff)):
            # direct x-great x-grand parent
            if i2 in my_parental_gen:
                relationship = (great-1) * 'great-' + (grand) * 'grand' + 'parent'
                return relationship
            # direct aunt/uncle
            # increase my parents by 1
            my_parental_gen = get_parents(family_graph, my_parental_gen, 1)
            if check_shared_parents(my_parental_gen, relative_parental_gen):
                relationship = (great) * 'great-' + 'aunt-uncle'
            # some kind of cousin
            else:
                for g in range(height_i1):
                    cousin += 1
                    my_parental_gen = get_parents(family_graph, my_parental_gen, 1)
                    relative_parental_gen = get_parents(family_graph, relative_parental_gen, 1)
                    if check_shared_parents(my_parental_gen, relative_parental_gen):
                        relationship = str(cousin) + '-cousin-' + str(removed) + '-removed'
                        return relationship
            grand = 1
            great += 1

    return relationship

def get_parents(family_graph, seeds, gen):
    if gen == 0:
        return family_graph[seeds[0]].parents
    for i in range(gen):
        ng_parents = []
        for parent in seeds:
            parents_parents = family_graph[parent].parents
            if parents_parents != []:
                for p in parents_parents:
                    ng_parents.append(p)
        seeds = ng_parents
    return ng_parents

def check_shared_parents(i1_parents, i2_parents):
    shared_parents = False
    for parent in i1_parents:
        if parent in i2_parents:
            return True
    for parent in i2_parents:
        if parent in i1_parents:
            return True
    return shared_parents

if __name__ == '__main__':
    main()
