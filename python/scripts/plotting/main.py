import readEncoding
import utils
import writeSegments
import plotOne
import plotTwo
import plotTwoB
import plotR

segmentLength = 10

encodedData = 'encoding.txt' #'All.wgs.svs.genotypes.encoded.txt1
segCounts = './segcounts/' + str(segmentLength) + 'segcounts.txt'
freqCounts = './freqcounts/' + str(segmentLength) + 'freqcounts.txt'
pltOneName = './png/' + str(segmentLength) + 'freqloglog.png'
pltTwoName = './png/segcounts.png'
pltTwoBName = './png/numelements.png'

freqRtxt = './segcounts/Big10segcounts.txt'
plotRName = './png/' + str(segmentLength) + 'segcounts.png'

segCountsDir = '/Users/kristen/PycharmProjects/precision-medicine/segcounts/'

def main():
    print('Reading encoded data.\n')
    all_lines = readEncoding.file_to_lists(encodedData)
    print('Making segments.\n')
    segmented_lines = readEncoding.segment_lines(all_lines, segmentLength)
    print('Counting unique segments.\n')
    unique_segments_per_region = utils.count_unique_segemnts_in_region(segmented_lines)
    unique_segment_count = utils.count_num_unique_segments(unique_segments_per_region)
    x = 'break'
    print('Writing \n')
    writeSegments.write_unique_segment_counts(unique_segment_count, segCounts)
    # unique_segments = utils.find_unique_segments(segmented_lines)
    # unique_seg_counts = utils.count_unique_segments(unique_segments, segmented_lines)
    # print('Writing segment counts.\n')
    # writeSegments.write_seg_count(unique_seg_counts, segCounts)
    # print('Writing frequency counts.\n')
    # writeSegments.write_freq(segCounts, freqCounts)
    #print('Plotting.\n')
    #plotR.plot_frequencey(segCounts, plotRName)
    # # plotOne.plot_frequencey(freqCounts, pltOneName)
    # plotTwo.plot_numElements(segCountsDir, pltTwoName)
    # plotTwoB.plot_numElements(encodedData, pltTwoBName)

    # print(str(len(unique_segments)) + " unique segments with a segment size of " + str(segmentLength) + ".")

if __name__ == '__main__':
    main()

