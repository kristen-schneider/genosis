#include <iostream>
#include <cctype>
#include <algorithm>
#include <fstream>
#include <string>

#include "readEncoding.h"

using namespace std;

float* read_encodings(string encodingtxt, int numSamples, int numVariants){
	
	// path to encoding file
	ifstream inFile;
	
	// open encoded data (.txt file)
        inFile.open(encodingtxt);
	if ( !inFile.is_open() ) {
                cout << "Failed to open: " << encodingtxt << endl;
        }


	// array to store all samples
	// arr[s][v]
	float* fArr = new float [numSamples * numVariants];

	string s;
	float f;
	int v  = 0;
	// read file, line by line
	if(inFile.is_open()){
		while(true){
			inFile >> s;
			if(inFile.eof()){ break; }
			f = stof(s);
			fArr[v] = f;
			v++;
		}
   	}
	cout << "total variants: " << v << endl;
	inFile.close();
	
	//// print out array
	//for (int i = 0; i < 3; i ++){
	//	for (int j = 0; j < 9; j++){
	//		cout << sampleArr[i][j];
	//	}
	//	cout << endl;
	//}

	return fArr;
}

// read queries file
float* read_queries(string queriestxt, int numVariants, int numQueries){
        // path to queries file
        ifstream inFile;

        // open queries data (.txt file)
        inFile.open(queriestxt);
        if ( !inFile.is_open() ) {
                cout << "Failed to open: " << queriestxt << endl;
        }

	// to store all queries
        float* queriesArr = new float[numVariants * numQueries];
	int Q = 0;
	if(inFile.is_open()){
		string line;
                while(getline(inFile, line)){
			int segLength = line.length();
			string s;
                        float f;

                        // convert string line to float array
                        for (int c = 0; c < segLength; c++){
                                s = line[c];
                                f = stof(s);
                                queriesArr[Q * segLength + c] = f;
                        }
			Q++;
		}
	}

        inFile.close();

        return queriesArr;

}
