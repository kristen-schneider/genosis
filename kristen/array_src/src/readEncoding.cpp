#include <iostream>
#include <fstream>
#include <string>

#include "readEncoding.h"

using namespace std;

float* read_test(int numSamples, int numVariants){
	
	// path to encoding file
        string inFileString = "/home/sdp/precision-medicine/encoding.txt";
	ifstream inFile;
	
	// open encoded data (.txt file)
        inFile.open(inFileString);
	if ( !inFile.is_open() ) {
                cout << "Failed to open: " << inFileString << endl;
        }


	// array to store all samples
	// arr[s][v]
	string sampleArr[numSamples][numVariants];
	float* fArr = new float [numSamples * numVariants];

	string s;
	float f;
	// read file, line by line
	if(inFile.is_open()){

		for (int i = 0; i < (numSamples * numVariants); i ++){
			inFile >> s;
			f = stof(s);
			fArr[i] = f;
		}
	//	for(int s = 0; s < numSamples; s++){
        //    		for(int v = 0; v < numVariants; v++){
	//			inFile >> sampleArr[s][v];
	//		}
        //	}
   	}

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
