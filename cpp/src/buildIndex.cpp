#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "buildIndex.h"


/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;

//faiss::IndexHNSWFlat build_faiss_index(string encodedTXT, int num_variants, int num_samples){
faiss::IndexFlatL2 build_faiss_index(string encodedTXT, int num_variants, int num_samples, char delim){
	// setup for FAISS
	faiss::IndexFlatL2 index(num_variants);
	//faiss::IndexHNSWFlat index(num_variants, 64);

	// open encoded file
	ifstream eFile;
	eFile.open(encodedTXT);
	if (!eFile.is_open()){
		cout << "Failed to open: " << encodedTXT << endl;
	}else{
		string line;
		while (getline(eFile, line)){
			string s;
			float f;

			// convert line to float array
			int i = 0;
			float* sample_vector = new float[num_variants];
			size_t start;
			size_t end = 0;
			if (delim == '\0'){
			 	for (int i = 0; i < num_variants; i ++){
					s = line[i];
					f = stof(s);
					sample_vector[i] = f;
				}
			}
			else{
				while ((start = line.find_first_not_of(delim, end)) != std::string::npos){
					end = line.find(delim, start);
					f = stof(line.substr(start, end - start));
					sample_vector[i] = f;
					i ++;
			
				}
			}
			/*
			for (int i = 0; i < num_variants; i ++){
				cout << sample_vector[i] << " ";
			}
			cout << endl;
			*/
			index.add(1, sample_vector);
			delete[] sample_vector;
		}
	}
	cout << "...added " << index.ntotal << " samples to the index." << endl;
	eFile.seekg(0);
	eFile.close();
	eFile.clear();
	return index;
}


//faiss::IndexHNSWFlat build_faiss_index_segments(string encodedFile, int start, int lengthSegment, int numSamples){
////faiss::IndexFlatL2 build_faiss_index_segments(string encodedFile, int start, int lengthSegment, int numSamples){
//	cout << "INDEX_HNSW_FLAT" << endl;
//	
//	// setup for FAISS
//	faiss::IndexHNSWFlat IndexHNSWFlat(lengthSegment, 64);
//	//faiss::IndexFlatL2 index(lengthSegment);
//	
//	//if (index.is_trained == 1){cout << "...index is trained." << endl;}
//	//else{cerr << "...INDEX IS NOT TRAINED." << endl;}
//	
//	// ifstream to encoded file
//        ifstream inFile;
//	// open encoded file
//        inFile.open(encodedFile);
//        if ( !inFile.is_open() ) {
//                cout << "Failed to open: " << encodedFile << endl;
//        }
//
//	// read encoded file line by line
//	string line;
//	int lineCount = 0;
//	if(inFile.is_open()){
//                while(getline(inFile, line)){
//			string s;
//			float f;
//			// convert string line to float array
//			float* singleVector = new float[lengthSegment];
//			int i = 0;
//			for (int c = start; c < start+lengthSegment; c++){
//				s = line[c];
//				f = stof(s);
//				singleVector[i] = f;	
//				i++;
//			}
//			
//			IndexHNSWFlat.add(1, singleVector);	
//			delete[] singleVector;
//			lineCount++;
//		}
//
//	}
//	cout << "...added " << IndexHNSWFlat.ntotal << " vectors to index." << endl;
//	// closed encoded file
//	inFile.seekg(0);
//	inFile.close();
//	inFile.clear();
//	return IndexHNSWFlat;
//}
