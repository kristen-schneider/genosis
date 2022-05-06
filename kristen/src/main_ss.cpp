#include<iostream>
#include<string>
#include<vector>

//#include "read_encodings.h"
//#include "faiss_pm.h"
#include "utils.h"

using namespace std;

//g++ main_ss.cpp read_encodings.cpp faiss_pm.cpp -I/home/sdp/miniconda3/envs/py38/include/ -L/home/sdp/miniconda3/envs/py38/lib/ -lfaiss -o test ss
int main(void) {

	cout << "Starting program: Part 2." << endl;

	cout << "...reading encoded file..." << endl;
	int numSamples = 3; // to make space for an array of strings which holds a string encoding for each sample
	int numVariants = 9;
	string cohort_arr[numSamples]; // to store all samples as a list of strings
	vector<float> vecOfFloats;
	//vecOfFloats = read_encoded_data(numSamples);
	
	cout << "...starting similarity search..." << endl;
	//float cohort_float_arr[numSamples];
	//faiss_flat(numSamples, numVariants, cohort_arr, cohort_float_arr);

	cout << "End of program: Part 2." << endl;

	return 0;
}
