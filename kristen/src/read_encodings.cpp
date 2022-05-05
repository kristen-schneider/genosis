#include <iostream>
#include <fstream>
#include <string>

#include "read_encodings.h"

using namespace std;

void read_encoded_data(){

	// path to encoding file
        string inFileString = "../../encoding.txt";
	ifstream inFile;

        inFile.open(inFileString);
	if ( !inFile.is_open() ) {
                cout << "Failed to open: " << inFileString << endl;
        }
	
	string line;	// to store line from file
	if (inFile.is_open()) {
		while (getline (inFile, line)) {
			cout << line << std::endl;
		}
		inFile.close();
	}
}
