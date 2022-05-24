#include <iostream>
#include <fstream>
#include <string>

#include "chunkEncoding.h"

using namespace std;

float* split_encoding(string encodingtxt, int numSamples, int segmentLength){

	// path to encoding file
        ifstream inFile;
	// open encoded data (.txt file)
        inFile.open(encodingtxt);
        if ( !inFile.is_open() ) {
                cout << "Failed to open: " << encodingtxt << endl;
        }

	// array to store all samples
        float* fArr = new float [numSamples * segmentLength];
	
	string s;
	string segment;
        float f;
        // read file, line by line
        if(inFile.is_open()){

                for (int i = 0; i < (numSamples * segmentLength); i ++){
			inFile >> s;
			f = stof(s);
                        fArr[i] = f;
                }
	}
	inFile.close();	

	return fArr;
}
