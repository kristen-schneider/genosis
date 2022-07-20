## plots after talking to Ryan
def count_unique_segemnts_in_region(segments):
    numSegments = len(segments[0])
    unique = [[i] for i in segments[0]]

    for sample in range(len(segments)):
        curr_sample = segments[sample]
        for segment in range(numSegments):
            curr_segment = curr_sample[segment]
            if curr_segment not in unique[segment]:
                unique[segment].append(curr_segment)

    return unique

def count_num_unique_segments(unique_segments):
    unique_segment_count = dict()

    for u_s in unique_segments:
        num_unique_segments = len(u_s)
        try: unique_segment_count[num_unique_segments] += 1
        except KeyError: unique_segment_count[num_unique_segments] = 1
    return unique_segment_count


## Plots from Kristen
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