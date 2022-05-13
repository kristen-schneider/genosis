# Precision Medicine Project

## Release v 0.1: faiss search (basic) 
<br>
- install faiss with conda: <code>conda install -c pytorch faiss-cpu</code>
- compile with commandline: <code>g++ main_ss.cpp read_encodings.cpp faiss_pm.cpp -I /path/to/conda/include/ -I /../include/ -L /path/to/conda/lib/ -lfaiss -o faiss</code>
<br>
### NOTES:
- Variables that are currently hard coded: numSamples, numVariants, numQueries, xq, encodedtxt. These can be changed in main.<br>
- There is no error handling.
----------------------------------------------
<br><br>
## Kristen's Notes
<br>

### Separate Scripts...<br>
- Still having issues with .bashrc loading <code>LD_LIBRARY_PATH</code>. Must execute <code>export LD_LIBRARY_PATH=/usr/local/lib/:/usr/bin/:/home/sdp/miniconda3/envs/py38/lib/</code>.<br>
- FAISS installed with conda. Run FAISS part with with: <code>g++ main_ss.cpp read_encodings.cpp faiss_pm.cpp -I /home/sdp/miniconda3/envs/py38/include/ -L /home/sdp/miniconda3/envs/py38/lib/ -lfaiss -o test</code> --> still errors.<br>

### IDEALLY...<br>
I want to run <code>make</code> and my program would:<br>
1. Read VCF into Variant Major Format(VMF)<br>
2. Transpose VMF to Sample Major Format(SMF)<br>
3. Write out Sample Major Format to intermediate file<br>
4. Read intermediate file into arrays<br>
5. Use arrays as input to FAISS<br><br>


## PJ's Notes: <br>
<br>
mkdir build<br>
cd build<br>
cmake ..<br>
<br>
make
