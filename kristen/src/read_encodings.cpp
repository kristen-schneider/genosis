#include <iostream>
#include <iterator>
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
        string inFileString = "../../encoding_space.txt";
	ifstream inFile;
	// open encoded data (.txt file)
        inFile.open(inFileString);
	if ( !inFile.is_open() ) {
                cout << "Failed to open: " << inFileString << endl;
        }
	

	string line;	// to store line from file
	vector<vector<float>> vecVecOfFloats;	// to return at the end of function

	int varCount = 0;
	vector<float> vecOfFloats; // to store line as a vector of floats
	if (inFile.is_open()) {
		string s; 

		while(getline(inFile, s)){
		//string s = "What is the right way to split a string into a vector of strings";
			//vector<string> line_vector;
			stringstream ss(s);
			istream_iterator<string> begin(ss);
			istream_iterator<string> end;
			vector<string> vstrings(begin, end);
			//copy(vstrings.begin(), vstrings.end(), ostream_iterator<string>(cout, "\n"));
			
			// go through line vector and convert each entry to float
			for (int v = 0; v < vstrings.size(); v++){
				float f = stof(vstrings.at(v));
				vecOfFloats.push_back(f);
			}

			vecVecOfFloats.push_back(vecOfFloats);
			vecOfFloats.clear();
		}
		inFile.close();
	}

	return vecVecOfFloats;
}
