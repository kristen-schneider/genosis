#ifndef FAISS_UTILS_H
#define FAISS_UTILS_H

#endif //FAISS_UTILS_H

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <map>


using namespace std;

int count_num_samples(string in_file);
int count_length_input_vector(string in_file, char delim);
map<string, float*> make_ID_data_map(string data_file, char delim, int num_elements);
