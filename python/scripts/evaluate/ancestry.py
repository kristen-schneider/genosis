def label_ancestry(ancestry_file):
    super_ancestry = {}
    sub_ancestry = {}

    with open(ancestry_file, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            sample_ID = values[0]
            sub_pop = values[3]
            super_pop = values[5]

            super_ancestry[sample_ID] = super_pop
            sub_ancestry[sample_ID] = sub_pop

    file.close()
    return super_ancestry, sub_ancestry

def get_super_sub(ancestry_file):
    super_sub = {}

    with open(ancestry_file, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            sub_pop = values[3]
            super_pop = values[5]

            try:
                super_sub[super_pop].add(sub_pop)
            except KeyError:
                super_sub[super_pop] = set{}
                super_sub[super_pop].add(sub_pop)
    file.close()
    return super_sub

def get_ancestry_names(ancestry_file):
    super_labels = {}
    sub_labels = {}

    with open(ancestry_file, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            sub_ID = values[3]
            sub_name = values[4]
            super_ID = values[5]
            super_name = values[6]

            sub_labels[sub_ID] = sub_name
            super_labels[super_ID] = super_name
    file.close()
    return super_labels, sub_labels
