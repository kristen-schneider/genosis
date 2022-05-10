#include <iostream>
#include <fstream>
#include <vector>

#include "utils.h"
#include "slice.h"

using namespace std;

//g++ main_encode.cpp slice.cpp utils.cpp -lhts -o encode
int main(void){

	cout << "Start of encoding." << endl;
	
	cout << "...reading vcf and writing encodings to intermediate file..." << endl;
	sliceVCF();
	
	cout << "End of encoding." << endl;
	return 0;
}
