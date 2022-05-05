#include <iostream>
#include <fstream>
#include <vector>
#include "utils.h"
#include "slice.h"
#include "read_encodings.h"

using namespace std;

int main(void){

	cout << "Start of program." << endl;
	
	cout << "...reading vcf..." << endl;
	sliceVCF();
	
	cout << "...reading encoding..." << endl;
	read_encoded_data();

	cout << "End of program." << endl;
	return 0;
}
