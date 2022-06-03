# Precision Medicine Project

## Release v 0.1: faiss search (basic) 
<br>
- install faiss with conda: <code>conda install -c pytorch faiss-cpu</code><br>
- open <code>main.cpp</code> and make edits to variables in the top section.<br>
- compile with commandline: <code>g++ main_ss.cpp read_encodings.cpp faiss_pm.cpp -I /path/to/conda/include/ -I /../include/ -L /path/to/conda/lib/ -lfaiss -o faiss</code><br>

### NOTES: <br>
- Variables that are currently hard coded: numSamples, numVariants, numQueries, xq, encodedtxt. These can be changed in main.cpp.<br>
- There is no error handling.<br>
----------------------------------------------<br>

### Misc...<br>
- Still having issues with .bashrc loading <code>LD_LIBRARY_PATH</code>. Must execute <code>export LD_LIBRARY_PATH=/usr/local/lib/:/usr/bin/:/home/sdp/miniconda3/envs/py38/lib/</code>.<br>
- FAISS installed with conda. Run FAISS part with with: <code>g++ main_ss.cpp read_encodings.cpp faiss_pm.cpp -I /home/sdp/miniconda3/envs/py38/include/ -L /home/sdp/miniconda3/envs/py38/lib/ -lfaiss -o test</code> --> still errors.<br>
- to just encode: <code>g++ main.cpp readVCF.cpp utils.cpp -I ../include/ -lhts -o encode</code><br>

### WORKFLOW...<br>
! main.cpp scripts will need paths to proper input files. ! <br>
1. Read and encode VCF into Variant Major Format(VMF)<br>
<code>g++ main_encode.cpp readVCF.cpp utils.cpp -I ../include/ -lhts -o ../bin/encode</code><br>
2. Run FAISS.
<code>g++ main_faiss.cpp buildIndex.cpp searchIndex.cpp readEncoding.cpp -I ../include/ -I /home/sdp/miniconda3/include/ -L /home/sdp/miniconda3/lib/ -lfaiss -o ../bin/faiss</code><br>
3. Run comparison metrics. This includes running FAISS. (e.g. brute force euclidean distance vs FAISS)<br>
<code>g++ main_compare.cpp compare.cpp buildIndex.cpp searchIndex.cpp metrics.cpp readEncoding.cpp -I ../include/ -I /home/sdp/miniconda3/include/ -L /home/sdp/miniconda3/lib/ -lfaiss -o ../bin/compare</code><br>


## PJ's Notes: <br>
<br>
mkdir build<br>
cd build<br>
cmake ..<br>
<br>
make
