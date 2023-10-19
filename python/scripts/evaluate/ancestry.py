from collections import defaultdict

def get_ancestry_scores(samples, ancestry_labels, data_dir, ext):
    
    all_ancestry_scores = defaultdict(dict)
    
    for sample in samples:
        sample_file = data_dir + '/' + sample + ext + '.csv'
        query = sample
        q_ancestry = ancestry_labels[query]

        with open(sample_file, 'r') as file:
            # combine any sample matches across haps
            sample_haps = {}
            for line in file:
                if ":" not in line:
                    hap = line
                else:
                    values = line.strip().split()
                    match_ID = values[0].split(":")[0][:-2]
                    m_ancestry = ancestry_labels[match_ID]
                    score = float(values[1])
                    try:
                        sample_haps[match_ID].append(score)
                    except KeyError:
                        sample_haps[match_ID] = [score]
        # for all matches with this sample append and get score
        for m in sample_haps:
            curr_score = sum(sample_haps[m])
            relation = pair_relations[query][m]
            if relation == "sibling" and ext == '.pop':
                if curr_score > 100: print(query, m, sample_haps[m], curr_score)
            try:
                all_relations_scores[relation].append(curr_score)
            except KeyError:
                all_relations_scores[relation] = [curr_score]
    return all_relations_scores

def label_ancestry(ancestry_file):
    """
    returns dictionary with key: sample ID and value: population ID
    """
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
    """
    returns dictionary where key: superancestry and value: subancestry
    """
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
    """
    returns dictionary where key: ancestry ID and value: full ancestry name
    """
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
