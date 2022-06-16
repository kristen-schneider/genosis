#include <chrono>
#include <cstdlib>
#include <iostream>
#include <faiss/IndexFlat.h>
#include <fstream>
#include <sstream>

#include "buildIndex.h"
#include "searchIndex.h"

using namespace std;
using namespace std::chrono;
// 64-bit int
using idx_t = faiss::Index::idx_t;
template <class indexType>
indexType buildIndex(indexType index, string encodedFile, int start, int lengthSegment, int numSamples){
	if (index.is_trained == 1){cout << "...index is trained." << endl;}
        else{cerr << "...INDEX IS NOT TRAINED." << endl;}

        // ifstream to encoded file
        ifstream inFile;

        // open encoded file
        inFile.open(encodedFile);
        if ( !inFile.is_open() ) {
		cout << "Failed to open: " << encodedFile << endl;
        }

        // read encoded file line by line
        string line;
        int lineCount = 0;
        if(inFile.is_open()){
                while(getline(inFile, line)){
                        string s;
                        float f;
                        // convert string line to float array
                        float* singleVector = new float[lengthSegment];
                        int i = 0;
                        for (int c = start; c < start+lengthSegment; c++){
                                s = line[c];
                                f = stof(s);
                                singleVector[i] = f;
                                i++;
                        }
			index.add(1, singleVector);
                        delete[] singleVector;
                        lineCount++;
                }

        }

        cout << "...added " << index.ntotal << " vectors to index." << endl;
	inFile.close();
        inFile.seekg(0);
        inFile.clear();

        return index;
}

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(int argc, char* argv[]){

	// path to encoded file
	string encodingtxt = argv[1];		// file with encoded data
	string queriestxt = argv[2];		// file with query data

	int numVariants = atoi(argv[3]);	// number of variants (cols)
	int numSamples = atoi(argv[4]);		// number of samples (rows)
	int numQueries= atoi(argv[5]); 		// number of queries sumbitted
	int k = atoi(argv[6]); 			// num nearest neighbors
	int segmentLength = atoi(argv[7]); 	// length of one segment


	int numSegments = numVariants/segmentLength;
	cout << "NUM SEGMENTS (floor): " << numSegments << endl;


	// DONE. Start FAISS..
	cout << "---------" << endl;
	cout << "Starting similarity searching using FAISS..." << endl;
	
	cout << "\nENCODED FILE: " << encodingtxt << "..." << endl;
	//faiss::IndexFlatL2 index = build_faiss_index(encodingtxt, numVariants, numSamples);
	
	int start = 0;
	for (int i = 0; i < numSegments; i ++){
		auto startTime = high_resolution_clock::now();	
		cout << "\nSegment: " << start << "-" << start+segmentLength << endl;
		cout << "-Building index." << endl;
		faiss::IndexFlatL2 indexFL2(segmentLength);
        	//faiss::IndexFlatIP indexFIP(segmentLength);

        	buildIndex<faiss::IndexFlatL2>(indexFL2, encodingtxt, start, segmentLength, numSamples);
        	//buildIndex<faiss::IndexFlatIP>(indexFIP, encodingtxt, start, segmentLength, numSamples);
		//faiss::IndexFlatIP s_index = build_faiss_index_segments_IP(encodingtxt, start, segmentLength, numSamples);
		cout << "-Running similairty search." << endl;
		similarity_search(indexFL2, queriestxt, start, segmentLength, numVariants, numSamples, numQueries, k, to_string(start));
		start += segmentLength;
		
		auto stopTime = high_resolution_clock::now();
		auto durationTime = duration_cast<microseconds>(stopTime - startTime);
		cout << "TIME: " << durationTime.count() << endl;
	}
	if (numVariants % segmentLength != 0){
		int lastSegmentLength = numVariants - (numSegments * segmentLength);
		cout << "\nLAST SEG DIFF";// << endl;
		cout << "\nSegment: " << start << "-" << start+lastSegmentLength << endl;
                cout << "-Building index." << endl;
                faiss::IndexFlatL2 indexFL2(segmentLength);
        	//faiss::IndexFlatIP indexFIP(segmentLength);

        	buildIndex<faiss::IndexFlatL2>(indexFL2, encodingtxt, start, segmentLength, numSamples);
        	//buildIndex<faiss::IndexFlatIP>(indexFIP, encodingtxt, start, segmentLength, numSamples);

		//faiss::IndexFlatIP s_index = build_faiss_index_segments_IP(encodingtxt, start, lastSegmentLength, numSamples);
                cout << "-Running similairty search." << endl;
                similarity_search(indexFL2, queriestxt, start, lastSegmentLength, numVariants, numSamples, numQueries, k, to_string(start));

	}
	/*
	cout << "\n2.Running similairty search..." << endl;
	similarity_search(index, queriestxt, numVariants, numSamples, numQueries, k);	
	*/
	
	cout << "End of FAISS." << endl;
	cout << "---------" << endl;


	return 0;
}
