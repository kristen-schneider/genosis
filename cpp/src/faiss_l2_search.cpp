#include <cstdlib>
#include <faiss/IndexFlat.h>
#include <faiss/IndexHNSW.h>
#include <faiss/AutoTune.h>
#include <faiss/index_factory.h>
#include <faiss/index_io.h>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "faiss_l2_search.h"
#include "faiss_utils.h"


// 64-bit int
using idx_t = faiss::Index::idx_t;
using namespace std;

int main(int argc, char* argv[]){

        const char* index_file = argv[1];
	string database_IDs = argv[2];
	string query_IDs = argv[3];
	string query_encodings = argv[4];
	int k = stoi(argv[5]);
	string out_file = argv[6];

	ofstream out_faiss;
	out_faiss.open(out_file);

	char delim = ' ';
	int num_q_samples = count_num_samples(query_IDs);
	cout << "Number of queries: " << num_q_samples << endl;
	int vector_size = count_length_input_vector(query_encodings, delim);
	cout << "Length of query vector: " << vector_size << endl;

	vector<string> database_ID_vector = make_query_ID_vector(database_IDs);
	vector<string> query_ID_vector = make_query_ID_vector(query_IDs);

	cout << "Searching index..." << endl;
	faiss::Index* index = faiss::read_index(index_file);
	idx_t* I = new idx_t[k * num_q_samples];
	float* D = new float[k * num_q_samples];

	float* queries = make_queries_arr(query_encodings, query_IDs, delim, num_q_samples, vector_size);
	///vector<string> database_ID_vector = make_query_ID_vector(database_IDs);
	//vector<string> query_ID_vector = make_query_ID_vector(query_IDs);
	index->search(num_q_samples, queries, k, D, I);

	for (int i = 0; i < num_q_samples; i++){
		out_faiss << "QUERY: " << query_ID_vector.at(i) << endl;
		for (int j = 0; j < k; j++){
			out_faiss << database_ID_vector.at(I[i * k + j]) << " " << I[i * k + j] << "\t" << sqrt(D[i * k + j]) << endl;
		}
		out_faiss << endl;
	}	
	out_faiss.close();
}



/*
 * makes an array of all queries
 */
float* make_queries_arr(string query_data, string query_IDs, char delim, int num_q_samples, int num_elements){
	float* queries = new float[num_q_samples * num_elements];

	// make map of sample ID to sample encoding (embedding)
        map<string, float*> ID_query_map = make_ID_data_map(query_data, delim, num_elements);

        // read through query file (list of samples in database)
        // and pull those vectors from the map to add to index
        ifstream q_file_stream;
        q_file_stream.open(query_IDs);
        if (!q_file_stream.is_open()){
                cout << "Failed to open: " << query_IDs << endl;
                exit(1);
        }
        string line;
	int sample_i = 0;
        while (getline(q_file_stream, line)){
                // pull sample vector from map
                float* sample_vector = ID_query_map[line];
		// add all elements to query array
		for (int q = 0; q < num_elements; q++){
			queries[sample_i * num_elements + q] = sample_vector[q];
		}
		sample_i ++;
        }
        q_file_stream.close();
	cout << "finished making query array" << endl;
	return queries;
}

/*
 * Make a list of query IDs to print
 * for output.
 */
vector<string> make_query_ID_vector(string query_IDs){
	vector<string> query_ID_vector;

	// read through query ID file (list of samples in query)
        // and add each ID to a vector of string
        ifstream q_file_stream;
        q_file_stream.open(query_IDs);
        if (!q_file_stream.is_open()){
                cout << "Failed to open: " << query_IDs << endl;
                exit(1);
        }
        string line;
        while (getline(q_file_stream, line)){
		query_ID_vector.push_back(line);
	}

	return query_ID_vector;
}
