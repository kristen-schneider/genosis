def find_unique_segments(segments):
    unique = []

    for s in segments:
        for s_i in s:
            if s_i not in unique:
                unique.append(s_i)
    return unique


def count_unique_segments(unique_segments, all_segments):
    seg_counts = {k: 0 for k in unique_segments}
    for s in all_segments:
        for s_i in s:
            seg_counts[s_i] += 1

    return seg_counts