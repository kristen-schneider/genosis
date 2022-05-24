#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <faiss/IndexFlat.h>

#include "faiss_pm.h"
#include "bruteForce.h"

using namespace std;
/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */


// 64-bit int
using idx_t = faiss::Index::idx_t;


int ss(float* database, float* queries, int numSamples, int numVariants, int numQueries){ 
	//int d = numVariants;      // dimension
	//int nb = numSamples; // database size
	//int nq = numQueries;  // nb of queries


	faiss::IndexFlatL2 index(numVariants); // call constructor
	printf("is_trained = %s\n", index.is_trained ? "true" : "false");
	index.add(numSamples, database); // add vectors to the index
	printf("ntotal = %zd\n", index.ntotal);

	// debugging
//	for (int i = 0; i < ( numVariants * numSamples); i++){
//		printf("%f ", database[i]);
//	}
//	cout << endl << typeid(database[0]).name() << endl;
//	for (int i = 0; i < ( numVariants * numQueries); i++){
//		printf("%f ", queries[i]);
//	}
//	cout << endl << typeid(queries[0]).name() << endl;
	
	int k = 4; // number of nearest neightbors to return
	
	// sanity check
	{
		idx_t* I = new idx_t[k * numQueries];
                float* D = new float[k * numQueries];

		index.search(numQueries, database, k, D, I);

		// print results
                cout << "I=\n" << endl;
                for (int i = 0; i < numQueries; i++){
                        for (int j = 0; j < k; j++){
                                cout << "    " << I[i * k + j] << " ";
			}
                        cout << endl;
                }
		cout << "D=\n" << endl;
                for (int i = 0; i < numQueries; i++){
                        for (int j = 0; j < k; j++){
                                cout << "    " << D[i * k + j] << " ";
			}
                        cout << endl;
        	}

		// test accuracy
		FAISS_vs_BF(database, queries, numSamples, numVariants, numQueries, I, D, k);

		delete[] I;
        	delete[] D;

	}
	
//	// FAISS on database with queries
//	{
//		idx_t* I = new idx_t[k * numQueries];
//		float* D = new float[k * numQueries];
//
//		index.search(numQueries, queries, k, D, I);
//
//		// print results
//		cout << "RESULTS FOR I_SS:" << endl;
//		for (int i = 0; i < numQueries; i++){
//			cout << "  Query " << i << ": ";
//			for (int j = 0; j < k; j++){
//				cout << I[i * k + j] << "\t";
//			}
//			cout << endl;
//		}
//		cout << "RESULTS FOR D_SS:" << endl;
//		for (int i = 0; i < numQueries; i++){
//			cout << "  Query " << i << ": ";
//			for (int j = 0; j < k; j++){
//				cout << D[i * k + j] << "\t";
//			}
//			cout << endl;
//		}
//		
//		delete[] I;
//		delete[] D;
//	}
	
	delete[] database;
	delete[] queries;

	return 0;
}


/*
 *  to compare the distance given by FAISS to distance computed by Brute Force for k nearest neighbors
 * */
float FAISS_vs_BF(float* database, float* queries, int numSamples, int numVariants, int numQueries, idx_t* I, float* D, int k){
	
	cout << "Testing against brute force." << endl;

	float diff = 0;
	float eucDist = 0;

	// for k nearest neighbors returned
	int start = 0;
	int end = numVariants;
	for (int i = 0; i < numQueries; i++){
		float* query_slice = arrSlice(queries, start, end);

		for(int j = 0; j < k; j++){
		
			int i_index = I[i * k + j];
			const float* db_slice = arrSlice(database, i_index, i_index + numVariants);
			
			eucDist = euclidean_distance(db_slice, query_slice, numVariants);
			cout << eucDist << endl;
			
		}
		cout << endl;

		//float* db_slice = arrSlice(database, start, end);
		//float* query_slice = arrSlice(queries, start, end);	

		start = end;
		end += numVariants;

	}

	return diff;

}

