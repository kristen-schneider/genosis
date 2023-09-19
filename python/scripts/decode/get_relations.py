import argparse
from collections import defaultdict
import graph_helpers as gh
import tree_helpers as th

def parse_args():
    parser = argparse.ArgumentParser(description='Label pedigree relations between individuals')
    parser.add_argument('--ped', type=str, help='PED file')
    parser.add_argument('--root_p', type=str, help='Paternal Root Node')
    parser.add_argument('--root_m', type=str, help='Maternal Root Node')
    return parser.parse_args()

def main():
    args = parse_args()
    ped_file = args.ped
    relations_file = ped_file.replace('.fam', '.relations')
    root_p = args.root_p
    root_m = args.root_m

    family_tree = th.build_family_tree(ped_file)
    family_graph = gh.build_family_graph(ped_file)

    # gh.plot_tree(family_tree)
    family_members = family_tree.keys()
    relations_labels = defaultdict(dict)
    unrelated = False

    for i1 in family_members:
        for i2 in family_members:

            if i1 == "kristen" and i2 == "chloe":
                x = 0

            distance = gh.search_family_graph(family_graph, family_members, i1, i2, root_p, root_m)
            height_1 = th.get_node_height(i1, family_tree, family_members, root_p, root_m)
            height_2 = th.get_node_height(i2, family_tree, family_members, root_p, root_m)
            # for all nxn pairs of individuals, find their relationship
            unrelated = (th.detect_unrelated(i1, family_tree, root_p, root_m) or
                         th.detect_unrelated(i2, family_tree, root_p, root_m))
            if unrelated:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
                unrelated = False
            else:
                relation, reverse_relation = label_relations(i1, i2,
                                                             distance, height_1, height_2,
                                                             root_p, root_m)
            try:
                relations_labels[i1][i2] = relation
            except KeyError:
                relations_labels[i1] = {i2: relation}
            try:
                relations_labels[i2][i1] = reverse_relation
            except KeyError:
                relations_labels[i2] = {i1: reverse_relation}


    # write to file
    with open(relations_file, 'w') as file:
        for i1 in family_members:
            for i2 in family_members:
                file.write(i1 + ',' + i2 + ',' + relations_labels[i1][i2] + '\n')



def label_relations(i1, i2, distance, height_1, height_2, root_p, root_m):
    relation = ''
    reverse_relation = ''
    # can only be self
    if distance == 0 and height_1 == height_2 and i1 == i2:
        relation = 'self'
        reverse_relation = 'self'
    # can be a parent or child
    elif distance == 1:
        if height_1 > height_2:
            relation = 'child'
            reverse_relation = 'parent'
        elif height_2 > height_1:
            relation = 'parent'
            reverse_relation = 'child'
        else:
            relation = 'unrelated'
            reverse_relation = 'unrelated'
    # can be a grandparent or grandchild or sibling
    elif distance == 2:
        if (i1 == root_p or i1 == root_m) and (i2 == root_m or i2 == root_p):
            relation = 'unrelated'
            reverse_relation = 'unrelated'
        elif height_1 > height_2:
            relation = 'grandchild'
            reverse_relation = 'grandparent'
        elif height_2 > height_1:
            relation = 'grandparent'
            reverse_relation = 'grandchild'
        else:
            relation = 'sibling'
            reverse_relation = 'sibling'
    # can be a great-grandparent or
    # great-grandchild or
    # aunt/uncle or
    # niece/nephew
    elif distance == 3:
        if height_1 > height_2:
            if height_1 - height_2 == 1:
                relation = 'niece/nephew'
                reverse_relation = 'aunt/uncle'
            elif height_1 - height_2 == 2:
                relation = 'great-grandchild'
                reverse_relation = 'great-grandparent'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
        elif height_2 > height_1:
            if height_2 - height_1 == 1:
                relation = 'aunt/uncle'
                reverse_relation = 'niece/nephew'
            elif height_2 - height_1 == 3:
                relation = 'great-grandparent'
                reverse_relation = 'great-grandchild'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
        else:
            relation = 'unrelated'
            reverse_relation = 'unrelated'
    # can be a great-great-grandparent or
    # great-great-grandchild or
    # great-aunt/uncle or
    # great-niece/nephew or
    # 1st cousin
    elif distance == 4:
        if height_1 > height_2:
            if height_1 - height_2 == 3:
                relation = 'great-great-grandparent'
                reverse_relation = 'great-great-grandchild'
            elif height_1 - height_2 == 2:
                relation = 'great-aunt/uncle'
                reverse_relation = 'great-niece/nephew'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
        elif height_2 > height_1:
            if height_2 - height_1 == 1:
                relation = 'great-great-grandchild'
                reverse_relation = 'great-great-grandparent'
            elif height_2 - height_1 == 3:
                relation = 'great-niece/nephew/uncle'
                reverse_relation = 'great-aunt/uncle'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
        elif(height_1 == height_2):
            relation = '1st cousin'
            reverse_relation = '1st cousin'
        else:
            relation = 'unrelated'
            reverse_relation = 'unrelated'

    elif distance == 5:
        if height_1 > height_2:
            if height_1 - height_2 == 4:
                relation = 'great-great-great-grandparent'
                reverse_relation = 'great-great-great-grandchild'
            elif height_1 - height_2 == 3:
                relation = 'great-great-aunt/uncle'
                reverse_relation = 'great-great-niece/nephew'
            elif height_1 - height_2 == 1:
                relation = '1st-cousin-once-removed'
                reverse_relation = '1st-cousin-once-removed'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
        elif height_2 > height_1:
            if height_2 - height_1 == 4:
                relation = 'great-great-great-grandchild'
                reverse_relation = 'great-great-great-grandparent'
            elif height_2 - height_1 == 3:
                relation = 'great-great-niece/nephew'
                reverse_relation = 'great-great-aunt/uncle'
            elif height_2 - height_1 == 1:
                relation = '1st-cousin-once-removed'
                reverse_relation = '1st-cousin-once-removed'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'

        else:
            relation = 'unrelated'
            reverse_relation = 'unrelated'

    elif distance == -1:
        relation = 'unrelated'
        reverse_relation = 'unrelated'
    else:
        relation = 'undetermined'
    return relation, reverse_relation



if __name__ == '__main__':
    main()
