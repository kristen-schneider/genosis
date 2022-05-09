#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "read_encodings.h"


using namespace std;

// need to convert this array workflow to vectors
// returns vector of vector of floats for input to ss
vector<vector<float>> read_encoded_data(int numSamples, int numVariants){

	// path to encoding file
        string inFileString = "../../encoding.txt";
	ifstream inFile;
	// open encoded data (.txt file)
        inFile.open(inFileString);
	if ( !inFile.is_open() ) {
                cout << "Failed to open: " << inFileString << endl;
        }
	

	string line;	// to store line from file
	vector<vector<float>> vecVecOfFloats;	// to return at the end of function

	int varCount = 0;
	if (inFile.is_open()) {
		string tmp; 
		vector<float> vecOfFloats; // to store line as a vector of floats
		vector<string> words;

		while(getline(inFile, tmp, ',')){
			if(varCount < numVariants){
				float f = stof(tmp);
                        	vecOfFloats.push_back(f);
				varCount += 1;
			}else{
				cout << "line" << endl;
				vecVecOfFloats.push_back(vecOfFloats);
				varCount = 0;
			}
		}
		inFile.close();
	}
	cout << vecVecOfFloats.size() << endl;
	return vecVecOfFloats;
}
