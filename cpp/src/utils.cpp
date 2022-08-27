#include <iostream>
#include <fstream>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <string>
#include <vector>
#include <map>

//#include "utils.h"

using namespace std;

// transpose a vector of vectors of ints
vector<vector<int>> transpose(vector<vector<int>> &vmf){
        /*
         * Takes a variant major format
         * vector of vector of ints
         * and transposes it to
         * sample major format.
         */
        cout << "Transposing VMF to SMF..." << endl;
        // throw error if size is bad
        if (vmf.size() == 0) {
                cerr << "Error reading variant major format." << endl;
        }

        // transpose data
        vector<vector<int>> smf_vec(vmf[0].size(), vector<int>());
        for (size_t i = 0; i < vmf.size(); i++) {
                for (size_t j = 0; j < vmf[i].size(); j++) {
                        smf_vec[j].push_back(vmf[i][j]);
                }
        }
    return smf_vec;
}

/*
 * get a list of all sample IDs
 */
vector<string> get_sample_IDs(string sample_IDs_file){

	vector<string> sample_IDs_vec;
	
	// open file and check success
	ifstream sample_IDs_file_stream;
	sample_IDs_file_stream.open(sample_IDs_file);
	if ( !sample_IDs_file_stream.is_open() ){
		cout << "FAILED TO OPEN: " << sample_IDs_file << endl;
		exit(1);
	}
	string line;
	while(getline(sample_IDs_file_stream, line)){
		sample_IDs_vec.push_back(line);
	}

	return sample_IDs_vec;
}

/*
 * returns number of samples in VCF file
 */
int get_num_samples(bcf_hdr_t *vcf_header){

        int num_samples = -1;
        num_samples = bcf_hdr_nsamples(vcf_header);
        return num_samples;
}

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
	return num_variants;
}
