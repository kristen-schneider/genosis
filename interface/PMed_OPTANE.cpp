#include <fstream>
#include <vector>
#include <math.h>
#include <map>
#include <stdint.h>
#include <omp.h>
#include <chrono>
#include <stdio.h>
#include <vector>
#include <errno.h>
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include "htslib/synced_bcf_reader.h"
#include <libpmemobj++/container/concurrent_hash_map.hpp>
#include <libpmemobj++/container/string.hpp>
#include <libpmemobj++/container/basic_string.hpp>
#include <libpmemobj++/p.hpp>
#include <libpmemobj++/persistent_ptr.hpp>
#include <libpmemobj++/pool.hpp>
#include <libpmemobj++/allocator.hpp>
#include <libpmemobj++/make_persistent_atomic.hpp>

#include "slice.h"

typedef std::chrono::high_resolution_clock Clock;
using namespace pmem::obj; 
using vector_int_type = pmem::obj::vector<int>;
using hashmap_type = concurrent_hash_map<p <int>, pmem::obj::string >;
using vec_vec_type = pmem::obj::vector<vector_int_type>;
//using hashmap_type_iv = concurrent_hash_map<p <int>, vector_int_type >;
//using hashmap_type_iv2 = concurrent_hash_map<p <int>, vector_int_type >;

struct root_pmem {
	persistent_ptr<hashmap_type> pptr;
    persistent_ptr<vec_vec_type> pptr_gt;
    persistent_ptr<vec_vec_type> pptr_gtT;

};
struct timespec begin1, end1, begin0, end0, begin11, end11, begin00, end00;
timespec diff(timespec start, timespec end )
{
    timespec temp;

    if ((end.tv_nsec-start.tv_nsec)<0)
    {
            temp.tv_sec = (end.tv_sec-start.tv_sec-1);
            temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
    }
    else 
    {
            temp.tv_sec = (end.tv_sec-start.tv_sec);
            temp.tv_nsec = (end.tv_nsec-start.tv_nsec);
    }
    return temp;
}
std::vector<std::vector<int> > transpose( std::vector<std::vector<int> > &b)
{
    if (b.size() == 0)
        //return ;
        std::cerr << "Error reaching to db." << std::endl;

    std::vector<std::vector<int>> trans_vec(b[1].size(), std::vector<int>());
    std::cout << b.size() << " \t::: \t"<< b[1].size() << std::endl;
    for (size_t i = 0; i < b.size(); i++)
    {
        for (size_t j = 0; j < b[i].size(); j++)
        {
            trans_vec[j].push_back(b[i][j]);
        }
    }
    std::cout << trans_vec.size() << " \t::: \t"<< trans_vec[1].size() << std::endl;
    return trans_vec;
}
//pmempool create obj --layout="node_hash_map" --size 20G /mnt/Optane/poolvsf_V01
//UNIVERSE OF DATA
// 0|0 -> 0
// 0|1, 1|0, 0|2, 2|0, 0|3, 3|0, 1|2, 2|1, 1|3, 3|1 ->1
// 1|1, 2|2, 3|3 ->2
// .|., 0|., .|0, 1|., .|1, 2|., .|2, 3|., .|3 ->3

int main(){

	// path to encoding file
        std::ifstream inFile;
        inFile.open("encoding.txt");
	
	std::string line;	// to store line from file
	if (inFile.is_open()) {
		while (getline (inFile, line)) {
			std::cout << line << std::endl;
		}
		inFile.close();
	}
	return 0;
}

//int main(){
//	std::cout << "Starting to run program..." << std::endl;
//	std::string s = "...testing sliceVCF...";
//	sliceVCF(s);
//}

int writeEncoding(void){
	std::map<std::string,int > datU = {{"0|0", 0}, {"0|1", 1}, {"1|0", 1}, {"0|2", 1}, {"2|0", 1},
     	{"0|3", 1}, {"3|0", 1}, {"1|2", 1}, {"2|1", 1}, {"1|3", 1}, {"3|1", 1}, 
    	{"1|1", 2},{"2|2", 2},{"3|3", 2},{".|.", 3},{"0|.", 3},{".|0", 3},{"1|.", 3},
     	{".|1", 3},{"2|.", 3},{".|2", 3},{"3|.", 3},{".|3", 3}};

	// counters
	int nn   = 0;  // total number of records in file
	int nsnp = 0;  // number of SNP records in file
	int nhq  = 0;  // number of SNPs for the single sample passing filters
	int nseq = 0;  // number of sequences
	
	// path to out file
	std::ofstream outFile;
	outFile.open("encoding.txt");

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
            printf("Loading codes...\n");
	
	std::vector<std::vector<int>> tempVecVec; // vector of vectors to storee all genotype encodings
	// read the VCF records, one by one
	while(bcf_read(test_vcf, test_header, test_record) == 0){
		bcf_unpack(test_record, BCF_UN_ALL);
		bcf_unpack(test_record, BCF_UN_INFO);
		
		std::string chrom = bcf_hdr_id2name(test_header, test_record->rid); // #CHROM
		int pos = (unsigned long)test_record->pos; // POS
		std::string tid =   test_record->d.id ; // ID
                std::string ref = test_record->d.allele[0]; // REF
                std::string alt = test_record->d.allele[1]; // ALT
                std::double_t qual =   test_record->qual; // QUAL
		
		// CHROM thru QUAL as one vector of strings.
		std::vector<std::string> strvec = 
		{chrom,(std::to_string(pos)),tid,ref,alt,(std::to_string(qual))};
		// CHROM thru QUAL as one string
		std::string s;
		for (const auto &piece : strvec) s += "\t" + piece;
		std::string strs = chrom +"\t" + 
			(std::to_string(pos)) + "\t" + 
			tid + "\t" + 
			ref + "\t" + 
			alt+ "\t" + 
			(std::to_string(qual));


		// genotypes
		std::string kg;
		std::vector<int> tempVec;
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
				kg = (std::to_string(genotype))+"|"+(std::to_string(genotype_s));
			}

			// one record is tempVec
			tempVec.push_back(datU[kg]);
		}
		// all records are tempVecVec
		tempVecVec.push_back(tempVec);
		for(int i=0; i < tempVec.size(); i++) {
			outFile << tempVec.at(i) << " ";
			//std::cout << tempVec.at(i) << " ";
		}
		tempVec.clear(); 
		//std::cout << "------------------" << std::endl;
	
		outFile << "\n";
		//std::cout << "\n";
	} // end of reading records

	// end of script
	return 0;
}


//int main(){
//    pool<root_pmem> pop;
//    //std::map<std::string,int > datU = {{"0|0", 0}, {"0|1", 1}, {"1|0", 1}, {"0|2", 1}, {"2|0", 1},
//    // {"0|3", 1}, {"3|0", 1}, {"1|2", 1}, {"2|1", 1}, {"1|3", 1}, {"3|1", 1}, 
//    //{"1|1", 2},{"2|2", 2},{"3|3", 2},{".|.", 3},{"0|.", 3},{".|0", 3},{"1|.", 3},
//    // {".|1", 3},{"2|.", 3},{".|2", 3},{"3|.", 3},{".|3", 3}};
//
//    std::map<std::string,int > datU = {{"0|0", 5}, {"0|1", 5}, {"1|0", 5}, {"0|2", 5}, {"2|0", 5},
//     {"0|3", 5}, {"3|0", 5}, {"1|2", 5}, {"2|1", 1}, {"1|3", 1}, {"3|1", 1}, 
//     {"1|1", 5},{"2|2", 5},{"3|3", 5},{".|.", 5},{"0|.", 5},{".|0", 5},{"1|.", 5},
//     {".|1", 5},{"2|.", 5},{".|2", 5},{"3|.", 5},{".|3", 5}};
//    // counters
//    int nn   = 0;  // total number of records in file
//    int nsnp = 0;  // number of SNP records in file
//    int nhq  = 0;  // number of SNPs for the single sample passing filters
//    int nseq = 0;  // number of sequences
//    bool remove_hashmap = false;
//
//    try{
//        std::string path = "/mnt/Optane/poolvsf_V01";
//
//        
//        try {
//            pop = pool<root_pmem>::open(path, "node_hash_map");
//            
//        } catch (pmem::pool_error &e) {
//            std::cerr << e.what() << std::endl;
//            return -1;
//        }
//        auto &r = pop.root()->pptr;
//        auto &pgt = pop.root()->pptr_gt;
//        auto &pgtT = pop.root()->pptr_gtT;
//        if (r == nullptr) {
//            /* Logic when file was first opened. First, we have to
//             * allocate object of hashmap_type and attach it to the
//             * root object. */
//            pmem::obj::transaction::run(pop, [&] {
//                r = make_persistent<hashmap_type>();
//                pgt = make_persistent<vec_vec_type>();
//                pgtT = make_persistent<vec_vec_type>();
//            });
//            r->runtime_initialize();
//        }  else {
//            //cout << "here"<< endl;
//            /* Logic when hash_map already exists. After opening of
//             * the pool we have to call runtime_initialize()
//             * function in order to recalculate mask and check for
//             * consistency. */
//
//            r->runtime_initialize();
//            //nd->runtime_initialize();
//
//            /* Defragment the whole pool at the beginning. */
//            try {
//                r->defragment();
//                //nd->defragment();
//            } catch (const pmem::defrag_error &e) {
//                std::cerr << "Defragmentation exception: "
//                      << e.what() << std::endl;
//                throw;
//            } 
//        }
//
//    
//    
//        auto &mapmt = *r;
//        auto &mapgt = *pgt;
//        auto &mapgtT = *pgtT;
//        uint32_t N = mapmt.size();
//               //const char *VCFPath = "/home/sdp/precidion-medicine/data/ALL.wgs.svs.genotypes.vcf";
//               const char *VCFPath = "/home/sdp/precision-medicine/data/short.vcf";
//        htsFile *test_vcf = bcf_open(VCFPath, "r");
//        if ( !test_vcf ) {
//            printf("Failed to open: %s\n", VCFPath);
//            }
//        // returning a bcf_hdr_t struct 
//        bcf_hdr_t *test_header = bcf_hdr_read(test_vcf);
//        fprintf(stderr, "File %s contains %i samples\n", "short.vcf", bcf_hdr_nsamples(test_header));
//        if(test_header == NULL) {throw std::runtime_error("Unable to read header.");}
//        const char **seqnames = NULL;
//        seqnames = bcf_hdr_seqnames(test_header, &nseq);
//        
//	// initialize and allocate bcf1_t object
//        bcf1_t *test_record = bcf_init();
//         if (test_record == NULL) {
//            fprintf(stderr, "ERROR: record is empty\n");
//        }
//        if (N==0){
//		printf("N = 0\n");
//            // genotype data for each call
//            // genotype arrays are twice as large as
//            // the other arrays as there are two values for each sample
//            int ngt_arr = 0;
//            int *gt     = NULL;
//            int ngt     = 0;
//            printf("Loading codes... ");
//            fflush(stdout);
//            //#pragma omp parallel for
//            //for (int t=0;t<parsvec.size();t++){
//            //    mapp.insert_or_assign(t, parsvec[t]);
//                //vecArray.push_back(parsvec[t].data());
//
//            //}
//            /*
//            for (int i = 0; i < nseq; i++) {
//                // bcf_hdr_id2name is another way to get the name of a sequence
//                fprintf(stderr, "[%2i] %s (bcf_hdr_id2name -> %s)\n", i, seqnames[i],
//                        bcf_hdr_id2name(test_header, i));
//                std::cout << "------------------" << std::endl;
//                std::cout << nseq << std::endl;
//            }
//            */
//            int nsample = bcf_hdr_nsamples(test_header);
//            printf("%03" PRId32 "\n",nsample);
//            std::vector<std::vector<int>> tempVecVec ;
//            while(bcf_read(test_vcf, test_header, test_record) == 0){
//                
//                nn++;
//                //std::cout<< nn << std::endl;
//                bcf_unpack(test_record, BCF_UN_ALL);
//                bcf_unpack(test_record, BCF_UN_INFO);
//                //read CHROM
//                
//                std::string chrom = bcf_hdr_id2name(test_header, test_record->rid);
//                int pos = (unsigned long)test_record->pos ;
//                std::string tid =   test_record->d.id ;
//                std::string ref = test_record->d.allele[0];
//                std::string alt = test_record->d.allele[1];
//                std::double_t qual =   test_record->qual ;
//                std::vector<std::string> strvec = {chrom,(std::to_string(pos)),tid,ref,alt,(std::to_string(qual))};
//                std::string s;
//                for (const auto &piece : strvec) s += "\t" + piece;
//                std::string strs = chrom +"\t" +(std::to_string(pos))+"\t" +tid+"\t" +ref+"\t" +alt+"\t" +(std::to_string(qual));
//
//                mapmt.insert_or_assign(nn, strs);
//                //printf("nn is %d\n", nn);
//                std::string kg;
//                std::vector<int> tempVec;
//                ngt =  bcf_get_genotypes(test_header, test_record,  &gt, &ngt_arr);
//                int ngts = ngt/nsample;
//                for (int i=0; i<nsample; i++){
//                    //for (int j=0; j<ngts; j++){    
//                    int genotype_s;
//                    int genotype = bcf_gt_allele(gt[i*ngts+0]);    
//                        
//                     
//                    //}
//                    if (genotype==-1)
//                        kg = ".|.";
//                    else{
//                        genotype_s = bcf_gt_allele(gt[i*ngts+1]);
//                        kg = (std::to_string(genotype))+"|"+(std::to_string(genotype_s));
//                    }
//                    //std::cout<< datU[kg] << "\t"; 
//
//                    tempVec.push_back(datU[kg]);
//                }  
//                
//                mapgt.emplace_back(tempVec); 
//                tempVecVec.push_back(tempVec);
//                
//                  
//                tempVec.clear();      
//                //std::cout << "------------------" << std::endl; 
//                
//                ////check if is snp
//                //printf("ttt %d\n", bcf_is_snp(test_record));
//                
//            }
//            std::vector<std::vector<int>> TransposeTempVecVec = transpose(tempVecVec);
//            std::cout << "----|--------|------" << std::endl;
//            for (auto Tv:TransposeTempVecVec)
//                mapgtT.emplace_back(Tv);
//            free(gt);
//            free(seqnames);
//            bcf_hdr_destroy(test_header);
//            bcf_close(test_vcf);
//            bcf_destroy(test_record);            //
//            printf("Loaded\n");
//        }
//        else{
//            printf("reading codes from Pmem... \n");
//
//            for (int i = 1; i < 2; ++i) {
//                std::string tmpstr;
//                try{
//                    hashmap_type::accessor acc;
//                    bool mtres = mapmt.find(acc, i);
//                    
//                    if (mtres) {
//                        assert(acc->first == i);
//                        std::cout << i << "\t\t";
//                        std::cout << acc->second.c_str() <<std::endl;
//
//                    }
//                }
//                catch(int e){
//                    std::cout<< "ERROR" <<e <<std::endl;
//                }                    
//            }            
//            for (int i = 1; i < 2; ++i) {
//                //std::vector<std::string> tmpvect;
//                printf("SECOND for loop\n");
//		    try{
//                    std::cout << i << "\t\t";
//                    for (auto v:mapgt[i]) {                        
//                        std::cout << v  ;
//                    }
//                    std::cout << "" <<std::endl;
//                }
//                catch(int e){
//                    std::cout<< "ERROR" <<e <<std::endl;
//                }                    
//            }            
//            for (int j = 1; j < 2; ++j) {
//                //std::vector<std::string> tmpvect;
//		printf("THIRD for loop\n");
//		    try{
//                    std::cout << j << "\t\t";
//                    for (auto vT:mapgtT[j]) {
//                        std::cout << vT  ;
//                    }
//                    std::cout << "" <<std::endl;
//                }
//                catch(int e){
//                    std::cout<< "ERROR" <<e <<std::endl;
//                }                    
//            }  
//        }
//
//        try {
//            /* Defragment the whole pool at the end. */
//            mapmt.defragment();
//        } catch (const pmem::defrag_error &e) {
//            std::cerr << "Defragmentation exception: " << e.what()
//                << std::endl;
//            throw;
//        } catch (const pmem::lock_error &e) {
//            std::cerr << "Defragmentation exception: " << e.what()
//                << std::endl;
//            throw;
//        } catch (const std::range_error &e) {
//            std::cerr << "Defragmentation exception: " << e.what()
//                << std::endl;
//            throw;
//        } catch (const std::runtime_error &e) {
//            std::cerr << "Defragmentation exception: " << e.what()
//                << std::endl;
//            throw;
//        }
//
//        if (remove_hashmap) {
//            /* Firstly, erase remaining items in the map. This
//            * function is not thread-safe, hence the function is
//            * being called only after thread execution has
//            * completed. */
//            try {
//                mapmt.clear();
//                mapgt.clear();
//                mapgtT.clear();
//            } catch (const pmem::transaction_out_of_memory &e) {
//                std::cerr << "Clear exception: " << e.what()
//                    << std::endl;
//                throw;
//            } catch (const pmem::transaction_free_error &e) {
//                std::cerr << "Clear exception: " << e.what()
//                    << std::endl;
//                throw;
//            }
//
//            /* If hash map is to be removed, free_data() method
//            * should be called first. Otherwise, if deallocating
//            * internal hash map metadata in a destructor fails
//            * program might terminate. */
//            mapmt.free_data();
//            mapgt.free_data();
//            mapgtT.free_data();
//            /* map.clear() // WRONG
//            * After free_data() hash map cannot be used anymore! */
//
//            transaction::run(pop, [&] {
//                delete_persistent<hashmap_type>(r);
//                delete_persistent<vec_vec_type>(pgt);
//                delete_persistent<vec_vec_type>(pgtT);
//                r = nullptr;
//                pgt = nullptr;
//                pgtT = nullptr;
//            });
//        }
//        } catch (const pmem::transaction_out_of_memory &e) {
//            std::cerr << "Exception occurred: ooot" << e.what() << std::endl;
//        } catch (std::exception &e) {
//            //std::cerr << "Exception occured: eee" << e.what() << std::endl;
//        try {
//            pop.close();
//        } catch (const std::logic_error &e) {
//            std::cerr << "Exception: " << e.what() << std::endl;
//        }
//        return -1;
//    }
//    return 0;      
//      
//}

