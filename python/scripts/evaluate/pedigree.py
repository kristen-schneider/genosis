def read_relations(relations_file):
    relations = {}
    with open(relations_file, 'r') as file:
        for line in file:
            values = line.strip().split(',')
            familyID = values[0]
            i1 = values[1]
            i2 = values[2]
            relation = values[3]
            try:
                relations[i1].update({i2: relation})
            except KeyError:
                relations[i1] = {i2: relation}
    return relations


def get_pedigree_scores(samples, pair_relations, data_dir, ext):

    all_relations_scores = {}

    for sample in samples:
        #print(sample)
        # open sample file
        sample_file = data_dir + '/' + sample + ext + '.csv'
        query = sample
        with open(sample_file, 'r') as file:
            # combine any sample amtches across haps
            sample_haps = {}
            for line in file:
                if ":" not in line:
                    hap = line
                else:
                    values = line.strip().split()
                    match_ID = values[0].split(":")[0][:-2]
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
