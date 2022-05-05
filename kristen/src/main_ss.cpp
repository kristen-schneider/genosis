#include<iostream>
#include<string>

#include "read_encodings.h"
#include "faiss_pm.h"

using namespace std;

//g++ main_ss.cpp read_encodings.cpp faiss_pm.cpp -I/home/sdp/miniconda3/envs/py38/include/ -L/home/sdp/miniconda3/envs/py38/lib/ -lfaiss -o test ss
int main(void) {

	cout << "Starting program: Part 2." << endl;

	cout << "...reading encoded file..." << endl;
	int numSamples = 3; // to make space for an array of strings which holds a string encoding for each sample
	string cohort_arr[numSamples]; // to store all samples as a list of strings
	read_encoded_data(numSamples, cohort_arr);
	
	cout << "...starting similarity search..." << endl;
	faiss_flat(numSamples, cohort_arr);

	cout << "End of program: Part 2." << endl;

	return 0;
}
