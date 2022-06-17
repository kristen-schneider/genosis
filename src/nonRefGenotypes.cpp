#include <iostream>
#include <fstream>
#include <map>

#include "readEncoding.h"
#include "nonRefGenotypes.h"

using namespace std;

map<int, int> count_non_reference_genotypes(string queriesTxt, string encodingTxt, int nQ, int nS, int nV){

	int nonReferenceCount = 0;
	
	float* queries = read_queries(queriesTxt, nV, nQ);
	map <int, int> samplecountMap = encoding_count(encodingTxt, queries, nS, nV);

	return samplecountMap;
}
//	float* encodings = read_queries(encodingTxt, nV, nS);
map<int, int> encoding_count(string encodingTxt, float* queryArr, int nS, int nV){
	// open encoded file
	ifstream eFile;
	eFile.open(encodingTxt);
	if ( !eFile.is_open() ) {
                cout << "Failed to open: " << encodingTxt << endl;
        }

	// to store all queries
        float* sampleArr = new float[nV];
        int S = 0;
	int count = 0;
	map<int, int> sampleCountMap;
        if(eFile.is_open()){
                string line;
                while(getline(eFile, line)){
                        int segLength = line.length();
                        string s;
                        float f;

                        // convert string line to float array
                        for (int c = 0; c < segLength; c++){
                                s = line[c];
                                f = stof(s);
                                sampleArr[c] = f;
                        }
			
			count = count_nonref_genotypes(S, nV, sampleArr, queryArr);
                        sampleCountMap[S] = count;
			S++;
                }
        }

        eFile.close();
        eFile.seekg(0);
        eFile.clear();
	return sampleCountMap;
}
	
int count_nonref_genotypes(int sample, int nV, float* S, float* Q){
	cout << "SAMPLE " << sample << ": ";
	int nonReferenceCount = 0;
	for (int v = 0; v < nV; v++){
		float qV = Q[v];
		float sV = S[v];
		// count shared non-reference genotype
		if ((sV != 0) && (qV != 0)){
			nonReferenceCount++;
		}
	}
	cout << nonReferenceCount << endl;

	return nonReferenceCount;
}
