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

map<string,vector<int>> make_biallelic_encoding_map(string encodingFile){
	 /*
         * read a file with genotype and its ecoding int
         * a map where the key is the genotype and the 
         * value is the bi-allelc encoding for that genotype
         */

        // a map whose key is the genotype and whose value is the encoding
        map<string, vector<int>> encoding_map;

	// a biallelic encoding is two bits
	vector<int> biallelic_encoding;

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
	    char delim = ' ';

            // split by delimiter
            while(getline(eFile, line)){
	    	vector<string> vec_line;
		split_line(line, delim, vec_line);
            	vector<int> g_encoding;
		
		string genotype = vec_line[0];
		g_encoding.push_back(stoi(vec_line[1]));
		g_encoding.push_back(stoi(vec_line[2]));

                // add genotype encoding to map
                //pair<string, vector<int>> p (genotype, g_encoding);
                //encoding_map.insert(p);
		encoding_map[genotype] = g_encoding;
            }

        }
        return encoding_map;
}

/*
 * Splits a string by some delimiter
 
void split_line(const string &s, char delim, vector<string> &elems){
    stringstream ss;
    ss.str(s);
    string item;
    while (getline(ss, item, delim)) {
            elems.push_back(item);
    }
}*/
