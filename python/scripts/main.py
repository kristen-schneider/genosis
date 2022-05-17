import readEncoding
import utils
import writeSegments
import plotOne
import plotTwo
import plotTwoB

encodedData = "All.wgs.svs.genotypes.encoded.txt"
segCounts = './segcounts/50segcounts.txt'
freqCounts = './freqcounts/50freqcounts.txt'
pltOneName = './png/50freqloglog.png'
pltTwoName = './png/segcounts.png'
pltTwoBName = './png/numelements.png'
segmentLength = 50
segCountsDir = '/Users/kristen/PycharmProjects/precision-medicine/segcounts/'

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
    # plotOne.plot_frequencey(freqCounts, pltOneName)
    plotTwo.plot_numElements(segCountsDir, pltTwoName)
    plotTwoB.plot_numElements(encodedData, pltTwoBName)

    # print(str(len(unique_segments)) + " unique segments with a segment size of " + str(segmentLength) + ".")

if __name__ == '__main__':
    main()

