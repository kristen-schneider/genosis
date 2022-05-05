#include <iostream>
#include <fstream>
#include <vector>
#include "utils.h"
#include "slice.h"
#include "read_encodings.h"
#include "faiss_pm.h"

using namespace std;

int main(void){

	cout << "Start of program." << endl;
	
	cout << "...reading vcf and writing encodings to intermediate file..." << endl;
	sliceVCF();
	
	cout << "...reading encoding..." << endl;
	read_encoded_data();

	cout << "...performing FAISS on encoded data..." << endl;
	faiss_flat();	

	cout << "End of program." << endl;
	return 0;
}
