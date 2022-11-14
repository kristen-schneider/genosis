#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

#include "map_encodings.h"
#include "utils.h"

using namespace std;

map<string,int> make_encoding_map(string encodingFile){
	/*
	 * read a file with genotype and its ecoding int
	 * a map where the key is the genotype and the 
	 * value is the encoding for that genotype
	 */

	// a map whose key is the genotype and whose value is the encoding
   	map<string, int> encoding_map;

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
   	        start = end + delim.length();
   	        g_encoding = stoi(line.substr(start));
   	        start = 0U;

   	        // add genotype encoding to map
   	        pair<string, int> p (genotype, g_encoding);
   	        encoding_map.insert(p);
   	    }

   	}
   	return encoding_map;
}
