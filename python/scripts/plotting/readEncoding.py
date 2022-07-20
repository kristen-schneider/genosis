def file_to_lists(f):
    o = open(f, 'r')
    all_lines = []

    for line in o:
        single_line = line.strip()
        all_lines.append(single_line)
    return all_lines

def segment_lines(all_lines, seg_length):
    all_line_segments = []

    for line in all_lines:
        line_segments = []
        segment = ''
        i = 0
        j = 0
        while i < len(line):
            if j < seg_length:
                if line[i] != ' ':
                    segment += line[i]
                    j += 1
            elif j == seg_length:
                line_segments.append(segment)
                if line[i] != ' ':
                    segment = line[i]
                    j = 1
                else:
                    segment = ''
                    j = 0
            i += 1
        if i == len(line):
            line_segments.append(segment)

        all_line_segments.append(line_segments)
    return all_line_segments




