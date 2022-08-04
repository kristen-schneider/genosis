#include <iostream>
#include <string>

#include "readVCF.h"
#include "utils.h"

using namespace std;

// code to read VCF and write to encoded file is commented out
// only code for reading encoded file and performing FAISS will be run
int main(int argc, char* argv[]){

        // path to encoded file
        const char *vcfFile = argv[1];        // input vcf file
        string encodedFile = argv[2];		// output encoded file 
	
	//const char *vcfFile = "/home/sdp/precision-medicine/data/vcf/ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf";
	//string encodedFile = "/home/sdp/precision-medicine/data/encodded/ALL.chr14.encoded";
	//const char *vcfFile = argv[1]; 		// VCF file
	//string encodedFile = argv[2]; 		// encoded file
	
	cout << "Start of encoding." << endl;

	cout << "Reading VCF file." << endl;
	sliceVCF(vcfFile, encodedFile);
		
	cout << "End of encoding." << endl;
	return 0;
}
