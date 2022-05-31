#include <iostream>
#include <cctype>
#include <fstream>
#include <string>

#include "readEncoding.h"

using namespace std;

// read queries file
float* read_queries(string queriestxt, int numVariants, int numQueries){
        // path to queries file
        ifstream qFile;

        // open queries data (.txt file)
        qFile.open(queriestxt);
        if ( !qFile.is_open() ) {
                cout << "Failed to open: " << queriestxt << endl;
        }

	// to store all queries
        float* queriesArr = new float[numVariants * numQueries];
	int Q = 0;
	if(qFile.is_open()){
		string line;
                while(getline(qFile, line)){
			int segLength = line.length();
			string s;
                        float f;

                        // convert string line to float array
                        for (int c = 0; c < segLength; c++){
                                s = line[c];
                                f = stof(s);
                                queriesArr[Q * segLength + c] = f;
                        }
			Q++;
		}
	}

        qFile.close();
        qFile.seekg(0);
        qFile.clear();
        return queriesArr;

}
