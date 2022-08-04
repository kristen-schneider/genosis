#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

#include "compare.h"
#include "metrics.h"
#include "readEncoding.h"

using namespace std;

void compare_main(string encodedFile, string queriesFile, int start, int lengthQuery, int numVariants, int numSamples, int numQueries, int numSegments, int metric){
	
	// get queries from file
	float* queries = read_queries_segment(queriesFile, start, numVariants, lengthQuery, numQueries);

	// for each query, do a separate computation	
	for(int q = 0; q < numQueries; q++){
		float* currQuery = new float[lengthQuery];
		for(int i = 0; i < (lengthQuery); i++){
			currQuery[i] = queries[q * lengthQuery + i];
		}
		cout << "Query " << q << endl;
		
		float* metricArr = compute_one_query(currQuery, encodedFile, start, lengthQuery, numVariants, numSamples, numQueries, metric);
			
		cout << "...computation complete." << endl;

        	cout << "...writing results." << endl;
        	ofstream outMetricFile;
        	outMetricFile.open("/home/sdp/precision-medicine/data/txt/" + to_string(metric) + ".Results." +to_string(start)+".txt");
        	for (int i = 0; i < numSamples; i++){
                	outMetricFile << i << "\t" << metricArr[i] << endl;
        	}
        	outMetricFile.close();

		//for (int s = 0; s < numSegments; s++){
		//	float* distArr = compute_one_query(currQuery, encodedFile, start, lengthQuery, numVariants, numSamples, numQueries);
		//}

		delete[] currQuery;
	}
}

float *compute_one_query(float* query, string encodedFile, int start, int segLength, int numVariants, int numSamples, int numQueries, int metric){
	// ifstream to encoded file
        ifstream inFile;
        // open encoded file
        inFile.open(encodedFile);
        if ( !inFile.is_open() ) {
                cout << "Failed to open: " << encodedFile << endl;
        }

	// store index and distance in my own hashmap
	float* metricArr = new float[numSamples];
        
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

			// compute metric by switch statment
			float singleMetric = -1;
			switch(metric) {
  				case 0:{
					// euclidean distance
					singleMetric = euclidean_distance(query, singleVector, segLength);
    					break;}
				case 1:{
					// count mismatches
					singleMetric = mismatch(query, singleVector, segLength);
					break;}
				case 2:{
					// count shared nonRef genotypes
					singleMetric = smart_mismatch(query, singleVector, segLength);
					break;}
				case 3:{
					// count shared nonRef genotypes
					singleMetric = sharedNRG(query, singleVector, segLength);
					break;}
				case 4:{
					// count shared nonRef genotypes with weights
					singleMetric = sharedNRGWeighted(query, singleVector, segLength);
				       	break;}
			}
			metricArr[lineCount] = singleMetric;
			lineCount ++;
		}
	}

	return metricArr;
}

