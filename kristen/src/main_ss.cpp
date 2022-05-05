#include<iostream>

#include "read_encodings.h"
#include "faiss_pm.h"

using namespace std;

//g++ main_ss.cpp read_encodings.cpp faiss_pm.cpp -I/home/sdp/miniconda3/envs/py38/include/ -L/home/sdp/miniconda3/envs/py38/lib/ -lfaiss -o test ss
int main(void) {

	cout << "Starting program: Part 2." << endl;

	cout << "...reading encoded file..." << endl;
	int numSamples = 3;
	read_encoded_data(numSamples);

	cout << "End of program: Part 2." << endl;

	return 0;
}
