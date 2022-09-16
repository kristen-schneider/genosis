vq = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0]
va = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
vb = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
vc = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
vd = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

def kristen(v1, v2, gaps_allowed):
    kd = 0
    running_score = 0
    num_consecutive_matches = 0
    gaps = 0
    all_scores = []

    size_vector = len(v1)
    for v in range(size_vector):
        if int(v1[v]) == 1 and int(v1[v]) == int(v2[v]):
            num_consecutive_matches += 1
            running_score += num_consecutive_matches
        else:
            gaps += 1
            if gaps > gaps_allowed:
                if running_score != 0:
                    all_scores.append(running_score)
                running_score = 0
                num_consecutive_matches = 0
                gaps = 0
            else:
                continue
    if running_score != 0:
        all_scores.append(running_score)
    #print(all_scores)
    kd = sum(all_scores)
    return kd

gaps_allowed = [0,1,2,5]
for g in gaps_allowed:
    print('GAPS ALLOWED: ', g)
    print('QQ', kristen(vq, vq, g))
    print('QA', kristen(vq, va, g))
    print('QB', kristen(vq, vb, g))
    print('QC', kristen(vq, vc, g))
    print('QD', kristen(vq, vd, g))

