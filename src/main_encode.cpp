#include <iostream>
#include <string>

#include "readVCF.h"
#include "utils.h"

using namespace std;

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(void){
	
	const char *vcfFile = argv[1]; 		// VCF file
	string encodedFile = argv[2]; 		// encoded file
	
	cout << "Start of encoding." << endl;

	cout << "Reading VCF file." << endl;
	sliceVCF(vcfFile, encodedFile);
		
	cout << "End of encoding." << endl;
	return 0;
}
