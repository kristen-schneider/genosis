from collections import defaultdict

def make_query_dict(queries_dict, faiss_file, num_segments, seg_idx):
    
    header = True

    f = open(faiss_file)
    for line in f:
        # start aggregation
        if 'QUERY' in line:
            header = False
            # start aggregation
            l = line.strip().split()
            query_sample_id = l[1]
            continue
        elif header:
            continue
        elif len(line.strip().split()) < 2:
            continue
        elif 'TIME' in line:
            break
        
        # for the following lines)
        l = line.strip().split()
        sample_ID = l[0]
        sample_idx = l[1]
        sample_dist = l[2]
    
        try:
            queries_dict[query_sample_id][sample_ID][seg_idx] = sample_dist
        except KeyError:
            queries_dict[query_sample_id].update({sample_ID: ['None'] * num_segments})
            queries_dict[query_sample_id][sample_ID][seg_idx] = sample_dist

    f.close()
