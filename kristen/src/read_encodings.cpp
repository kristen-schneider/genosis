#include <iostream>
#include <fstream>
#include <string>

#include "read_encodings.h"

using namespace std;

string* read_encoded_data(int numSamples){

	// path to encoding file
        string inFileString = "../../encoding.txt";
	ifstream inFile;

        inFile.open(inFileString);
	if ( !inFile.is_open() ) {
                cout << "Failed to open: " << inFileString << endl;
        }
	

	string line;	// to store line from file
	string cohort_arr[numSamples]; // to store all samples as a list of strings
	int s = 0;
	if (inFile.is_open()) {
		while (getline (inFile, line)) {
			// append line to the string array
			cohort_arr[s] = line;
			s += 1;
		}
		inFile.close();
	}
	return cohort_arr;
}
