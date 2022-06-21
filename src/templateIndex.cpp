#include <faiss/IndexFlat.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "templateIndex.h"


/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;



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
			/*cout << "adding vector: ";
                        for (int i = 0; i < segLength; i++){
                                cout << singleVector[i];
                        }*/

                        // add array to index
                        index.add(1, singleVector);
                        delete[] singleVector;
                        lineCount++;
                }

        }

	cout << "...added " << index.ntotal << " vectors to index." << endl;
        // closed encoded file
        inFile.close();
        inFile.seekg(0);
        inFile.clear();
		
	return index;
}


int main(){
	
	int start = 0;
	int lengthSegment = 3;
	int numSamples = 3;
	string encodedFile = "/home/sdp/precision-medicine/data/encoded/short.encoded.txt";

	faiss::IndexFlatL2 indexFL2(lengthSegment);
	faiss::IndexFlatIP indexFIP(lengthSegment);;

	buildIndex<faiss::IndexFlatL2>(indexFL2, encodedFile, start, lengthSegment, numSamples);
	buildIndex<faiss::IndexFlatIP>(indexFIP, encodedFile, start, lengthSegment, numSamples);

	return 0;
}

