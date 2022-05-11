#include <iostream>
#include <fstream>
#include <string>

#include "readEncoding.h"

using namespace std;

int read_test(int numSamples, int numVariants){
	
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
	string sampleArr[3][9];


	// read file, line by line
	int sampleCount = 0;
	if(inFile.is_open()){
		string myArray[9];

		for(int i = 0; i < 9; i++){
            		inFile >> sampleArr[sampleCount][i];
        	}
		sampleCount++;
		
		// print out array
		for (int i = 0; i < 3; i ++){
			for (int j = 0; j < 9; j++){
				cout << sampleArr[i][j];
			}
		}
   	}

	inFile.close();

	return 0;
}
