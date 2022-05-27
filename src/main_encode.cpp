#include <iostream>
#include <string>

#include "readVCF.h"
#include "utils.h"

using namespace std;

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(void){

	const char *vcfFile = "/home/sdp/precision-medicine/data/vcf/ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf";
	string encodedFile = "/home/sdp/precision-medicine/data/encoded/new.encoded.txt";
	cout << "Start of encoding." << endl;

	cout << "Reading VCF file." << endl;
	sliceVCF(vcfFile, encodedFile);

//	cout << "Reading Encoded file." << endl;
//	float* xb = read_test(encodingtxt, numSamples, numVariants);
//	cout << "Done Reading Encoded file." << endl;
//	
//
//	cout << endl << "Starting FAISS." << endl;
//	ss(xb, xq, numSamples, numVariants, numQueries);
//	cout << "End of FAISS." << endl;

	return 0;
}
