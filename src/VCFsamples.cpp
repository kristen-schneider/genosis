#include <iostream>
#include <fstream>
#include <math.h>
#include <map>
#include <string>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>

using namespace std;

int main(int argc, char **argv){
        const char *vcfFile = "/home/sdp/precision-medicine/data/vcf/ALL.chr14.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf";
// counters
        int nseq = 0;

        // open VCF/BCF file
        //    * use '-' for stdin
        //    * bcf_open will open bcf and vcf files
        //    * bcf_open is a macro that expands to hts_open
        //    * returns NULL when file could not be opened
        //    * by default also writes message to stderr if file could not be found
        htsFile * inf = bcf_open(argv[1], "r");
        if (inf == NULL) {
                return EXIT_FAILURE;
        }

        // read header
        bcf_hdr_t *hdr = bcf_hdr_read(inf);
        fprintf(stderr, "File %s contains %i samples\n", argv[1], bcf_hdr_nsamples(hdr));

        // report names of all the sequences in the VCF file
        const char **seqnames = NULL;
        // bcf_hdr_seqnames returns a newly allocated array of pointers to the seq names
        // caller has to deallocate the array, but not the seqnames themselves; the number
        // of sequences is stored in the int pointer passed in as the second argument.
        // The id in each record can be used to index into the array to obtain the sequence
        // name
        seqnames = bcf_hdr_seqnames(hdr, &nseq);
        fprintf(stderr, "Sequence names:\n");
        for (int i = 0; i < nseq; i++) {
                // bcf_hdr_id2name is another way to get the name of a sequence
                fprintf(stderr, "  [%2i] %s (bcf_hdr_id2name -> %s)\n", i, seqnames[i],
                       bcf_hdr_id2name(hdr, i));
        }


        // clean up memory
        if (seqnames != NULL)
                free(seqnames);
        bcf_hdr_destroy(hdr);
        bcf_close(inf);
        return 0;
}
