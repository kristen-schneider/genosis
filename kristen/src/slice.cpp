#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>

#include "slice.h"

using namespace std;

void sliceVCF(void){
        map<string,int > datU = {{"0|0", 0}, {"0|1", 1}, {"1|0", 1}, {"0|2", 1}, {"2|0", 1},
        {"0|3", 1}, {"3|0", 1}, {"1|2", 1}, {"2|1", 1}, {"1|3", 1}, {"3|1", 1},
        {"1|1", 2},{"2|2", 2},{"3|3", 2},{".|.", 3},{"0|.", 3},{".|0", 3},{"1|.", 3},
        {".|1", 3},{"2|.", 3},{".|2", 3},{"3|.", 3},{".|3", 3}};

        // counters
        int nn   = 0;  // total number of records in file
        int nsnp = 0;  // number of SNP records in file
        int nhq  = 0;  // number of SNPs for the single sample passing filters
        int nseq = 0;  // number of sequences

        // path to out file
        ofstream outFile;
        outFile.open("/home/sdp/precision-medicine/encoding.txt");

        // path to VCF file
        // short.vcf
        // ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf
        const char *VCFPath =
                "/home/sdp/precision-medicine/data/short.vcf";
	
	// open VCF file with htslib
        htsFile *test_vcf = bcf_open(VCFPath, "r");
      	if ( !test_vcf ) {
       		printf("Failed to open: %s\n", VCFPath);
       	}

	// returning a bcf_hdr_t struct
        bcf_hdr_t *test_header = bcf_hdr_read(test_vcf);
        int numSamples = 0;
        numSamples = bcf_hdr_nsamples(test_header); // counting number of samples

        fprintf(stderr, "File '%s' contains %i samples.\n", VCFPath, numSamples);
        if(test_header == NULL) {
                throw std::runtime_error("Unable to read header.");
        }
        const char **seqnames = NULL;
        seqnames = bcf_hdr_seqnames(test_header, &nseq); // getting sequence names

        // initialize and allocate bcf1_t object
        bcf1_t *test_record = bcf_init();
        if (test_record == NULL) {
                fprintf(stderr, "ERROR: record is empty\n");
        }

	// genotype data for each call
        // genotype arrays are twice as large as
        // the other arrays as there are two values for each sample
        int ngt_arr = 0;
        int *gt     = NULL;
        int ngt     = 0;
        printf("Encoding VCF to VMF...\n");

        vector<vector<int>> tempVecVec; // vector of vectors to storee all genotype encodings
}

