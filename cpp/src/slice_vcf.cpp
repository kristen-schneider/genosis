#include <iostream>
#include <fstream>
#include <string>

#include "slice_vcf.h"

using namespace std;

void write_vcf_header(string vcf_file, string out_file){
    /*
     * open VCF file and write header data to out_file
     */

    // open file and check success
    ifstream vcf_file_stream;
    vcf_file_stream.open(vcf_file);
    if (!vcf_file_stream.is_open()){
            cout << "FAILED TO OPEN: " << vcf_file << endl;
    }
    else{
        // open output VCF file for a segment
        ofstream out_file_stream;
        out_file_stream.open(out_file);

        string line;
        // read vcf file header and write just header to put
        while (getline (vcf_file_stream, line)){
            char char1 = line.at(0);
            if (char1 == '#'){
                out_file_stream << line << endl;
            }
            else{ break; }
        }
    }
}

int slice(string vcf_file, int segment_size, string base_name, string out_dir){
    /*
     * Open VCF file, ignore header,
     * write a segment of columns to out_file
     */

    // to return
    int segment_count = 0;

    // open vcf file and check success
    ifstream vcf_file_stream;
    vcf_file_stream.open(vcf_file);
    if (!vcf_file_stream.is_open()){
        cout << "FAILED TO OPEN: " << vcf_file << endl;
    }

    // read vcf file
    else{
        string line;
        int total_line_count = 0;
        int lines_in_segment = 0;

        // name segment out_file
        string out_vcf_segment_name = out_dir + base_name + \
                          ".seg." + to_string(segment_count) + \
                          ".vcf";
        // write vcf header
        write_vcf_header(vcf_file, out_vcf_segment_name);

        // open out file in append mode
        ofstream out_file_stream;
        out_file_stream.open(out_vcf_segment_name, ios_base::app);

        // read vcf file until segment length
        while (lines_in_segment < segment_size){
            if (vcf_file_stream.peek() != EOF){
                getline (vcf_file_stream, line);
            }
            else {
                lines_in_segment = segment_size;
                cout << "segment " << segment_count << " " << out_vcf_segment_name << endl;
                cout << "END OF FILE" << endl;
                //break;
            }

            // ingore header
            if (line.at(0) != '#'){
                //cout << total_line_count << " " \
                    << lines_in_segment << " " \
                    << out_vcf_segment_name << endl;
                out_file_stream << line << endl;
                // increment counters
                lines_in_segment ++;
                total_line_count ++;
            }
            // if segment length fulfilled
            if (lines_in_segment == segment_size){
                cout << "segment " << segment_count << " " << out_vcf_segment_name << endl;
                lines_in_segment = 0;
                out_file_stream.close();

                // open new file stream
                segment_count++;
                out_vcf_segment_name = out_dir + base_name + \
                                              ".seg." + to_string(segment_count) + \
                                              ".vcf";
                // write vcf header to new segment out_file
                write_vcf_header(vcf_file, out_vcf_segment_name);
                // open new segment out_file in append mode
                out_file_stream.open(out_vcf_segment_name, ios_base::app);

            }
        }
        out_file_stream.close();
        cout << "SEGMENT LENGTH: " << segment_size << endl;
        cout << "LINE COUNT: " << total_line_count << endl;
        cout << "NUM SEGMENTS: " << segment_count << endl;
        int num_variants = total_line_count;
	}
	return segment_count;
}
