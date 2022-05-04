#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
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
                throw runtime_error("Unable to read header.");
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

	cout << "testing" << endl;
	// read the VCF records, one by one
        while(bcf_read(test_vcf, test_header, test_record) == 0){
                bcf_unpack(test_record, BCF_UN_ALL);
                bcf_unpack(test_record, BCF_UN_INFO);

                string chrom = bcf_hdr_id2name(test_header, test_record->rid); // #CHROM
                int pos = (unsigned long)test_record->pos; // POS
                string tid =   test_record->d.id ; // ID
                string ref = test_record->d.allele[0]; // REF
                string alt = test_record->d.allele[1]; // ALT
                double_t qual =   test_record->qual; // QUAL

		// CHROM thru QUAL as one vector of strings.
		vector<string> strvec = 
		{chrom,(to_string(pos)),tid,ref,alt,(to_string(qual))};
		// CHROM thru QUAL as one string
		string s;
		for (const auto &piece : strvec) s += "\t" + piece;
		string strs = chrom +"\t" + 
			(to_string(pos)) + "\t" + 
			tid + "\t" + 
			ref + "\t" + 
			alt+ "\t" + 
			(to_string(qual));


		// genotypes
		string kg;
		vector<int> tempVec;
		ngt =  bcf_get_genotypes(test_header, test_record,  &gt, &ngt_arr);
		int ngts = ngt/numSamples;
		for (int i=0; i<numSamples; i++) {
			int genotype_s;
			int genotype = bcf_gt_allele(gt[i*ngts+0]);
			if (genotype==-1) {
				kg = ".|.";
			}
			else{
				genotype_s = bcf_gt_allele(gt[i*ngts+1]);
				kg = (to_string(genotype))+"|"+(to_string(genotype_s));
			}

			// one record is tempVec
			tempVec.push_back(datU[kg]);
		}
		// all records are tempVecVec
		tempVecVec.push_back(tempVec);
		//for(int i=0; i < tempVec.size(); i++) {
		//	outFile << tempVec.at(i) << " ";
		//	//cout << tempVec.at(i) << " ";
		//}
		tempVec.clear(); 
		//cout << "------------------" << endl;
	
		//outFile << "\n";
		//cout << "\n";

                cout << "record\n";
        } // end of reading records
}

