#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

#include "compare.h"
#include "metrics.h"
#include "readEncoding.h"
#include "sharedNRG.h"

using namespace std;

int compare_main(string encodedFile, string queriesFile, int start, int lengthQuery, int numVariants, int numSamples, int numQueries, int numSegments, string metric){
	
	// get queries from file
        //float* queries = read_queries(queriesFile, numVariants, numQueries);
	float* queries = read_queries_segment(queriesFile, start, numVariants, lengthQuery, numQueries);
	for(int q = 0; q < numQueries; q++){
		// one query at a time
		float* currQuery = new float[lengthQuery];
		for(int i = 0; i < (lengthQuery); i++){
			currQuery[i] = queries[q * lengthQuery + i];
		}
		cout << "Query " << q << endl;
		float* distArr = compute_one_query(currQuery, encodedFile, start, lengthQuery, numVariants, numSamples, numQueries, metric);

		//for (int s = 0; s < numSegments; s++){
		//	float* distArr = compute_one_query(currQuery, encodedFile, start, lengthQuery, numVariants, numSamples, numQueries);
		//}

		// sort distArr, keep indexes
		float* sortedDistArr = new float[numSamples];


		delete[] currQuery;
	}
	return 0;
}

float *compute_one_query(float* query, string encodedFile, int start, int segLength, int numVariants, int numSamples, int numQueries, string metric){
	// ifstream to encoded file
        ifstream inFile;
        // open encoded file
        inFile.open(encodedFile);
        if ( !inFile.is_open() ) {
                cout << "Failed to open: " << encodedFile << endl;
        }

	// store index and distance in my own hashmap
	float* distArr = new float[numSamples];
        
	// read encoded file line by line
        string line;
        int lineCount = 0;
        if(inFile.is_open()){
                while(getline(inFile, line)){
                        string s;
                        float f;

			// convert string line to float array
                        float* singleVector = new float[segLength];
                        int i = 0;
			for (int c = start; c < start+segLength; c++){
                                s = line[c];
                                f = stof(s);
                                singleVector[i] = f;
				i ++;
                        }
			float singleDistance = euclidean_distance(query, singleVector, segLength);
			distArr[lineCount] = singleDistance;
			lineCount ++;
		}
	}
	cout << "...brute force computations complete." << endl;

	// writing results
        cout << "...writing brute force results." << endl;
        ofstream outBruteForceFile;
        outBruteForceFile.open("/home/sdp/precision-medicine/data/txt/bruteforceResults." +to_string(start)+".txt");
        for (int i = 0; i < numSamples; i++){
                outBruteForceFile << i << "\t" << distArr[i] << endl;
        }
        outBruteForceFile.close();

	return distArr;
}

