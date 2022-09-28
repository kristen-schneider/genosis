vq = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0]
va = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
vb = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
vc = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
vd = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

def kristen_(v1, v2, gaps_allowed):
    kd = 0
    small_running_score = 0
    kd = 0
    num_consecutive_matches = 0
    gaps = 0
    
    size_vector = len(v1)
    for v in range(size_vector):
        if int(v1[v]) == 1 and int(v1[v]) == int(v2[v]):
            num_consecutive_matches += 1
            small_running_score += num_consecutive_matches
        else:
            gaps += 1
            if gaps > gaps_allowed:
                if small_running_score != 0:
                    kd += small_running_score
                small_running_score = 0
                num_consecutive_matches = 0
                gaps = 0
            else:
                continue
    if small_running_score != 0:
        kd += small_running_score
    
    # print(all_scores)
    return kd

def kristen(v1, v2, gaps_allowed):
    kd = 0
    small_running_score = 0
    num_consecutive_mismatches = 0
    gaps = 0

    size_vector = len(v1)
    for v in range(size_vector):
        #print(kd)
        if (int(v1[v]) != int(v2[v])) and (int(v1[v]) == 1 or int(v2[v]) == 1):        
            num_consecutive_mismatches += 1
            small_running_score += num_consecutive_mismatches
            
        else:
            gaps += 1
            if gaps > gaps_allowed:
                if small_running_score != 0:
                    kd += small_running_score
                small_running_score = 0
                num_consecutive_mismatches = 0
                gaps = 0
            else:
                continue
    if small_running_score != 0:
        kd += small_running_score

    # print(all_scores)
    return kd

# running program
gaps_allowed = [0,1,2,5]
for g in gaps_allowed:
    print('GAPS ALLOWED: ', g)
    print('QQ', kristen(vq, vq, g))
    print('QA', kristen(vq, va, g))
    print('QB', kristen(vq, vb, g))
    print('QC', kristen(vq, vc, g))
    print('QD', kristen(vq, vd, g))

