import argparse
from collections import defaultdict
import find_relations as fr

def parse_args():
    parser = argparse.ArgumentParser(description='Label pedigree relations between individuals')
    parser.add_argument('--ped', type=str, help='PED file')
    parser.add_argument('--out_dir', type=str, help='output directory')
    parser.add_argument('--roots', type=str, help='path to a list of root for all pedigrees')
    return parser.parse_args()

def main():
    args = parse_args()
    ped_file = args.ped
    out_dir = args.out_dir	
    relations_file = out_dir + 'samples.relations'
    roots = args.roots

    family_members = get_family_members(ped_file)
    family_roots = get_family_roots(roots)

    family_graph = fr.build_family_graph(ped_file)

    # gh.plot_tree(family_tree)
    relations_labels = defaultdict(dict)

    for i1 in family_members:
        roots_i1 = family_roots[family_members[i1]]
        if roots_i1 == []:
            continue
        for i2 in family_members:
            roots_i2 = family_roots[family_members[i2]]
            if roots_i2 == []:
                continue
            # check if same family
            if family_members[i1] == family_members[i2]:
                distance = fr.search_family_graph(family_graph, family_members, i1, i2, family_roots)
                # set height to 0 if i1 or i2 is root
                if i1 in roots_i1:
                    height_1 = 0
                else:
                    # height of i1 is i1s distance from root
                    height_1 = fr.search_family_graph(family_graph, family_members, i1, roots_i1[0], family_roots)
                if i2 in roots_i2:
                    height_2 = 0
                else:
                    # height of i2 is i2s distance from root
                    height_2 = fr.search_family_graph(family_graph, family_members, i2, roots_i2[0], family_roots)

                relation, reverse_relation = label_relations(i1, i2,
                                                             distance, height_1, height_2,
                                                             family_members, family_roots)
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'

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



def label_relations(i1, i2, distance, height_1, height_2, family_members, family_roots):

    roots_i1 = family_roots[family_members[i1]]
    roots_i2 = family_roots[family_members[i2]]

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
            relation = 'undetermined'
            reverse_relation = 'undetermined'
    # can be a grandparent or grandchild or sibling
    elif distance == 2:
        if (i1 in roots_i1) and (i2 in roots_i2):
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
            elif height_1 - height_2 == 3:
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
            relation = 'undetermined'
            reverse_relation = 'undetermined'
    # can be a great-great-grandparent or
    # great-great-grandchild or
    # great-aunt/uncle or
    # great-niece/nephew or
    # 1st cousin
    elif distance == 4:
        if height_1 > height_2:
            if height_1 - height_2 == 4:
                relation = 'great-great-grandchild'
                reverse_relation = 'great-great-grandparent'
            elif height_1 - height_2 == 2:
                relation = 'great-niece/nephew'
                reverse_relation = 'great-aunt/uncle'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
        elif height_2 > height_1:
            if height_2 - height_1 == 4:
                relation = 'great-great-grandparent'
                reverse_relation = 'great-great-grandchild'
            elif height_2 - height_1 == 2:
                relation = 'great-aunt/uncle'
                reverse_relation = 'great-niece/nephew'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
        elif(height_1 == height_2):
            relation = '1st cousin'
            reverse_relation = '1st cousin'
        else:
            relation = 'undetermined'
            reverse_relation = 'undetermined'

    elif distance == 5:
        if height_1 > height_2:
            if height_1 - height_2 == 5:
                relation = 'great-great-great-grandchild'
                reverse_relation = 'great-great-great-grandparent'
            elif height_1 - height_2 == 3:
                relation = 'great-great-niece/nephew'
                reverse_relation = 'great-great-auant/uncle'
            elif height_1 - height_2 == 1:
                relation = '1st-cousin-once-removed'
                reverse_relation = '1st-cousin-once-removed'
            else:
                relation = 'undetermined'
                reverse_relation = 'undetermined'
        elif height_2 > height_1:
            if height_2 - height_1 == 4:
                relation = 'great-great-great-grandparent'
                reverse_relation = 'great-great-great-grandchild'
            elif height_2 - height_1 == 3:
                relation = 'great-great-aunt/uncle'
                reverse_relation = 'great-great-niece/nephew'
            elif height_2 - height_1 == 1:
                relation = '1st-cousin-once-removed'
                reverse_relation = '1st-cousin-once-removed'
            else:
                relation = 'unrelated'
                reverse_relation = 'unrelated'
        else:
            relation = 'undetermined'
            reverse_relation = 'undetermined'
    elif distance == -1:
        relation = 'unrelated'
        reverse_relation = 'unrelated'
    elif distance == -2:
        relation = 'in-law-flag'
        reverse_relation = 'in-law-flag'
    else:
        relation = 'undetermined'
        reverse_relation = 'undetermined'

    return relation, reverse_relation

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
            roots = values[1:3]
            family_roots[family_ID] = roots
    return family_roots

if __name__ == '__main__':
    main()
