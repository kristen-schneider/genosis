#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <chrono>

#include "searchIndex.h"
#include "readEncoding.h"

/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std::chrono;
using namespace std;

//void search(const faiss::IndexHNSWFlat &index, int k, string queriesTXT,\
//	       	int num_queries, int num_variants){
void search(const faiss::IndexFlatL2 &index, int k, string queriesTXT,\
	       	int num_queries, int num_variants, const char delim){
	
	idx_t* I = new idx_t[k * num_queries];
	float* D = new float[k * num_queries];

	float* queries = read_queries(queriesTXT, num_variants, num_queries, delim); 
	
	auto start = high_resolution_clock::now();
	index.search(num_queries, queries, k, D, I);
	auto stop = high_resolution_clock::now();
	auto duration_search = duration_cast<microseconds>(stop - start);	
	cout << "TIME:search:" << duration_search.count() << endl;

	for (int i = 0; i < num_queries; i++){
		cout << "QUERY: " << i << endl;
		for (int j = 0; j < k; j++){
			cout << I[i * k + j] << "\t" << sqrt(D[i * k + j]) << endl;
		}
		cout << endl;
	}
	delete[] I;
	delete[] D;
}

void similarity_search(const faiss::IndexHNSWFlat &index, string queriesFile, int start, int lengthQuery, int numVariants, int numSamples, int numQueries, int k, string txtName){
//void similarity_search(const faiss::IndexFlatL2 &index, string queriesFile, int start, int lengthQuery, int numVariants, int numSamples, int numQueries, int k, string txtName){
	
	idx_t* I = new idx_t[k * numQueries];
	float* D = new float[k * numQueries];

	// get queries from file
	//float* queries = read_queries(queriesFile, lengthQuery, numQueries);
	cout << "...reading queries" << endl;
	float* queries = read_queries_segment(queriesFile, start, numVariants, lengthQuery, numQueries);
	/*
	cout << "Query: " << endl;
	for (int i = 0; i < lengthQuery; i++){
        	cout << queries[i];
        }
	cout << endl;
	*/

	//  index.search(nq, xq, k, D, I);
	
	index.search(numQueries, queries, k, D, I);
	//index.range_search(numQueries, queries, k, D, I, 50.0);
	cout << "...search complete." << endl;


	//// writing results
        //cout << "...writing index results." << endl;
	//ofstream outIndexFile;
	//outIndexFile.open("/home/sdp/precision-medicine/data/txt/indexResults.txt");
        for (int i = 0; i < numQueries; i++){
                for (int j = 0; j < k; j++){
                        cout << I[i * k + j] << "\t" << sqrt(D[i * k + j]) << endl;
                        //outIndexFile << I[i * k + j] << "\t" << sqrt(D[i * k + j]) << endl;
                }
                cout << endl;
                //outIndexFile << endl;
        }
	//outIndexFile.close();
	/*
        cout << "..writing distance results." << endl;
	ofstream outDistanceFile;
	outDistanceFile.open("/home/sdp/precision-medicine/data/txt/distanceReults.txt");
        for (int i = 0; i < numQueries; i++){
                for (int j = 0; j < k; j++){
                        outDistanceFile << "\t" << sqrt(D[i * k + j]) << " ";
                }
                outDistanceFile << endl;
        }
	outDistanceFile.close();
	*/
	
	delete [] I;
	delete [] D;
	cout << "...deleting I and D. " << "start: " << start << endl;
}

