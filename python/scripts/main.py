import readEncoding
import utils
import writeSegments
import plotOne

encodedData = "All.wgs.svs.genotypes.encoded.txt"
segCounts = '500segcounts.txt'
freqCounts = '500freqcounts.txt'
freqCountsOrd = '500freqcounts-ordered.txt'
pltName = '500freqloglog.png'
segmentLength = 200

def main():
    # print('Reading encoded data.\n')
    # all_lines = readEncoding.file_to_lists(encodedData)
    # print('Making segments.\n')
    # segmented_lines = readEncoding.segment_lines(all_lines, segmentLength)
    # print('Counting unique segments.\n')
    # unique_segments = utils.find_unique_segments(segmented_lines)
    # unique_seg_counts = utils.count_unique_segments(unique_segments, segmented_lines)
    # print('Writing segment counts.\n')
    # writeSegments.write_seg_count(unique_seg_counts, segCounts)
    # print('Writing frequency counts.\n')
    # writeSegments.write_freq(segCounts, freqCounts)
    print('Plotting.\n')
    plotOne.plot_frequencey(freqCountsOrd, pltName)
    # print(str(len(unique_segments)) + " unique segments with a segment size of " + str(segmentLength) + ".")

if __name__ == '__main__':
    main()

