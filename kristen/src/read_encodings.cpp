#include <iostream>
#include <fstream>

using namespace std;

void read_encoded_data(){

	// path to encoding file
        ifstream inFile;
        inFile.open("../encoding.txt");
	
	string line;	// to store line from file
	if (inFile.is_open()) {
		while (getline (inFile, line)) {
			cout << line << std::endl;
		}
		inFile.close();
	}
}
