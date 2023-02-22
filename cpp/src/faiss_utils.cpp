#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>

#include "faiss_utils.h"

using namespace std;

/*
 * Returns number of samples from a file with a 
 * list of samples.
 */
int count_num_samples(string in_file){
	int num_samples = 0;
	ifstream f_stream;
        f_stream.open(in_file);
        if ( !f_stream.is_open() ) {
                cout << "Failed to open: " << in_file << endl;
                exit(1);
        }
        // count number of lines = num_samples
        string line;
        while ( getline(f_stream, line) ){
		num_samples++;
	}
	f_stream.close();
	return num_samples;
}
/*
 * Returns number of entreis for a single
 * input vector  given an 
 * encoding or embedding file
 */
int count_length_input_vector(string in_file, char delim){
	int num_variants = 0;

	// open encoding or embeddding data
        ifstream f_stream;
        f_stream.open(in_file);
        if ( !f_stream.is_open() ) {
                cout << "Failed to open: " << in_file << endl;
                exit(1);
        }
        
	// count number of elements in first line = num_variants
        string line;
	size_t start;
        size_t end = 0;
	getline(f_stream, line);
	while ( (start = line.find_first_not_of(delim, end)) != std::string::npos){
		end = line.find(delim, start);
        	num_variants ++;
        }
	f_stream.close();
	// uncount sampleID
	return num_variants - 1;
}

map<string, float*> make_ID_data_map(string data_file, char delim, int num_elements){
        /*
         * Returns a map where key is sampleID and
         * value is a vector of floats for
         * the encoding (or embedding)
         */

        // open data_file
        ifstream data_file_stream;
        data_file_stream.open(data_file);
        if( !data_file_stream.is_open()){
                cout << "Failed to open: " << data_file << endl;
                exit(1);
        }

	// returned map
        map<string, float*> ID_data_map;

        cout << "mapping encoding vectors to samples" << endl;

        // read file
        string line;
        while (getline(data_file_stream, line)){
                size_t start;
                size_t end = 0;

                //vector<float> sample_vector;
                string sampleID;
		auto* sample_vector = new float[num_elements];
                int float_i = 0;
                while ((start = line.find_first_not_of(delim, end)) != std::string::npos){
                        end = line.find(delim, start);
                        if (sampleID.empty()){
                                sampleID = line.substr(start, end - start);
                        }else{
                                float f = stof(line.substr(start, end - start));
                                sample_vector[float_i] = f;
                                float_i++;
                        }
                }
		pair<string, float*> p (sampleID, sample_vector);
                ID_data_map.insert(p);
        }
        cout << "finished mapping encoding vectors" << endl;
        return ID_data_map;
}

