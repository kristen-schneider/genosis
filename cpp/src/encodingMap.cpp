#include <iostream>
#include <fstream>
#include <map>

using namespace std;

map<string,int> readEncodingFile(string encodingFile){
	// a map whose key is the genotype and whose value is the encoding
	map<string, int> encodingMap;
	
	// delimiter separates key and value
	string delim = " ";

	// open file with encoding specifications
	ifstream eFile;

	eFile.open(encodingFile);
	if ( !eFile.is_open() ) {
                cout << "Failed to open: " << encodingFile << endl;
        }
	
	// read file 
	if(eFile.is_open()){
		string line;
		string genotype;
		int g_encoding;
                
		// split by delimiter
		auto start = 0U;
		auto end = 0U;
		while(getline(eFile, line)){
			end = line.find(delim, start);
			genotype = line.substr(start, end);
			cout << genotype << endl;
			start = end + delim.length();
			end = line.find(delim, start);

		}
	}
	return encodingMap;
}

int validateEncodingCount(map<string,int> m){
	return m.size();
}

int main(void){
	string encodingFile = "/home/sdp/precision-medicine/encodings/trivial.txt";
	map<string, int> m = readEncodingFile(encodingFile);
	validateEncodingCount(m);
}
