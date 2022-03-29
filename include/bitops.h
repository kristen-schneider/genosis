/*
 * bitops.h
 *
 *  Created on: 2016-03-09
 *      Author: s2eghbal
 */

#ifndef BITOPS_H_
#define BITOPS_H_

#define popcntll __builtin_popcountll
#define popcnt __builtin_popcount
#define ctz __builtin_ctz

#include <stdio.h>
#include <math.h>
#include <iostream>
#include "types.h"
#include <bitset>

struct hammingR
{
	int r1;
	int r2;
};

const int lookup [] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
const double sqrt_lookup []= {0,1,1.4142,1.7321,2,2.2361,2.4495,2.6458,2.8284,3,3.1623,3.3166,3.4641,3.6056,3.7417,3.873,4,4.1231,4.2426,4.3589,4.4721,4.5826,4.6904,4.7958,4.899,5,5.099,5.1962,5.2915,5.3852,5.4772,5.5678,5.6569,5.7446,5.831,5.9161,6,6.0828,6.1644,6.245,6.3246,6.4031,6.4807,6.5574,6.6332,6.7082,6.7823,6.8557,6.9282,7,7.0711,7.1414,7.2111,7.2801,7.3485,7.4162,7.4833,7.5498,7.6158,7.6811,7.746,7.8102,7.874,7.9373,8};


inline int match(UINT8*P, UINT8*Q, int codelb) {
	switch(codelb) {
	case 4: // 32 bit
		return popcnt(*(UINT32*)P ^ *(UINT32*)Q);
	case 8: // 64 bit
		return popcntll(((UINT64*)P)[0] ^ ((UINT64*)Q)[0]);
	case 16: // 128 bit
        //std::cout<<popcntll(((UINT64*)P)[0] ^ ((UINT64*)Q)[0]) << "  "<< popcntll(((UINT64*)P)[1] ^ ((UINT64*)Q)[1])<<std::endl;
		return popcntll(((UINT64*)P)[0] ^ ((UINT64*)Q)[0]) \
				+ popcntll(((UINT64*)P)[1] ^ ((UINT64*)Q)[1]);
	case 32: // 256 bit
		return popcntll(((UINT64*)P)[0] ^ ((UINT64*)Q)[0]) \
				+ popcntll(((UINT64*)P)[1] ^ ((UINT64*)Q)[1]) \
				+ popcntll(((UINT64*)P)[2] ^ ((UINT64*)Q)[2]) \
				+ popcntll(((UINT64*)P)[3] ^ ((UINT64*)Q)[3]);
	case 64: // 512 bit
		return popcntll(((UINT64*)P)[0] ^ ((UINT64*)Q)[0]) \
				+ popcntll(((UINT64*)P)[1] ^ ((UINT64*)Q)[1]) \
				+ popcntll(((UINT64*)P)[2] ^ ((UINT64*)Q)[2]) \
				+ popcntll(((UINT64*)P)[3] ^ ((UINT64*)Q)[3]) \
				+ popcntll(((UINT64*)P)[4] ^ ((UINT64*)Q)[4]) \
				+ popcntll(((UINT64*)P)[5] ^ ((UINT64*)Q)[5]) \
				+ popcntll(((UINT64*)P)[6] ^ ((UINT64*)Q)[6]) \
				+ popcntll(((UINT64*)P)[7] ^ ((UINT64*)Q)[7]);
	default:
		int output = 0;
		for (int i=0; i<codelb; i++)
			output+= lookup[P[i] ^ Q[i]];
		return output;
	}
}
inline int popNOTAND(UINT8*P, UINT8*Q, int codelb)
{
	switch(codelb) {
	case 4: // 32 bit
		return popcnt(~(*(UINT32*)P) & *(UINT32*)Q);
	case 8: // 64 bit
	{
		return popcntll(~((UINT64*)P)[0] & ((UINT64*)Q)[0]);
	}
	case 16: // 128 bit
		return popcntll(~((UINT64*)P)[0] & ((UINT64*)Q)[0]) \
				+ popcntll(~((UINT64*)P)[1] & ((UINT64*)Q)[1]);
	case 32: // 256 bit
		return popcntll(~((UINT64*)P)[0] & ((UINT64*)Q)[0]) \
				+ popcntll(~((UINT64*)P)[1] & ((UINT64*)Q)[1]) \
				+ popcntll(~((UINT64*)P)[2] & ((UINT64*)Q)[2]) \
				+ popcntll(~((UINT64*)P)[3] & ((UINT64*)Q)[3]);
	case 64: // 512 bit
		return popcntll(~((UINT64*)P)[0] & ((UINT64*)Q)[0]) \
				+ popcntll(~((UINT64*)P)[1] & ((UINT64*)Q)[1]) \
				+ popcntll(~((UINT64*)P)[2] & ((UINT64*)Q)[2]) \
				+ popcntll(~((UINT64*)P)[3] & ((UINT64*)Q)[3]) \
				+ popcntll(~((UINT64*)P)[4] & ((UINT64*)Q)[4]) \
				+ popcntll(~((UINT64*)P)[5] & ((UINT64*)Q)[5]) \
				+ popcntll(~((UINT64*)P)[6] & ((UINT64*)Q)[6]) \
				+ popcntll(~((UINT64*)P)[7] & ((UINT64*)Q)[7]);
	default:
		int output = 0;
		for (int q=0; q<codelb; q++)
			output+= lookup[~P[q] & Q[q]];
		return output;
	}
}

inline double dotproduct(UINT8*P, UINT8*Q, int codelb)
{
	switch(codelb) {
	case 4: // 32 bit
		return popcnt((*(UINT32*)P) & *(UINT32*)Q);
	case 8: // 64 bi
		return popcntll(((UINT64*)P)[0] & ((UINT64*)Q)[0]);
	case 16: // 128 bit
		return popcntll(((UINT64*)P)[0] & ((UINT64*)Q)[0]) \
				+ popcntll(((UINT64*)P)[1] & ((UINT64*)Q)[1]);
	case 32: // 256 bit
		return popcntll(((UINT64*)P)[0] & ((UINT64*)Q)[0]) \
				+ popcntll(((UINT64*)P)[1] & ((UINT64*)Q)[1]) \
				+ popcntll(((UINT64*)P)[2] & ((UINT64*)Q)[2]) \
				+ popcntll(((UINT64*)P)[3] & ((UINT64*)Q)[3]);
	case 64: // 512 bit
		return popcntll(((UINT64*)P)[0] & ((UINT64*)Q)[0]) \
				+ popcntll(((UINT64*)P)[1] & ((UINT64*)Q)[1]) \
				+ popcntll(((UINT64*)P)[2] & ((UINT64*)Q)[2]) \
				+ popcntll(((UINT64*)P)[3] & ((UINT64*)Q)[3]) \
				+ popcntll(((UINT64*)P)[4] & ((UINT64*)Q)[4]) \
				+ popcntll(((UINT64*)P)[5] & ((UINT64*)Q)[5]) \
				+ popcntll(((UINT64*)P)[6] & ((UINT64*)Q)[6]) \
				+ popcntll(((UINT64*)P)[7] & ((UINT64*)Q)[7]);
	default:
		int output = 0;
		for (int q=0; q<codelb; q++)
			output+= lookup[~P[q] & Q[q]];
		return output;
	}
}

inline hammingR cosinematch(UINT8*P, UINT8*Q, int codelb) {
	hammingR R;
	R.r1 = popNOTAND(P,Q,codelb);
	//std::cout<<"R.r1 = "<<R.r1<<"\n";
	R.r2 = popNOTAND(Q,P,codelb);
	return R;
}
/* b <= 64 */
inline void split (UINT64 *chunks, UINT8 *code, int m, int mplus, int b) {
	UINT64 temp = 0x0;
	int nbits = 0;
	int nbyte = 0;
	UINT64 mask = b==64 ? 0xFFFFFFFFFFFFFFFFLLU : ((UINT64_1 << b) - UINT64_1);

	for (int i=0; i<m; i++) {
		while (nbits < b) {
			temp |= ((UINT64)code[nbyte++] << nbits);
			nbits += 8;
		}
		chunks[i] = temp & mask;
		temp = b==64 ? 0x0 : temp >> b;
		nbits -= b;
		if (i == mplus-1) {
			b--;		/* b <= 63 */
			mask = ((UINT64_1 << b) - UINT64_1);
		}
	}
}

/* generates the next binary code (in alphabetical order) with the
 * same number of ones as the input x. Taken from
 * http://www.geeksforgeeks.org/archives/10375
 */
inline UINT64 next_set_of_n_elements(UINT64 x) {
	UINT64 smallest, ripple, new_smallest;

	smallest     = x & -x;
	ripple       = x + smallest;
	new_smallest = x ^ ripple;
	new_smallest = new_smallest / smallest;
	new_smallest >>= 2;
	return ripple | new_smallest;
}

inline void print_code(UINT64 tmp, int b) {
	for (int j=(b-1); j>=0; j--) {
		printf("%llu", (long long int) tmp/(1 << j));
		tmp = tmp - (tmp/(1 << j)) * (1 << j);
	}
	printf("\n");
}

inline UINT64 choose(int n, int r) {
	UINT64 nchooser = 1;
	for (int k=0; k < r; k++) {
		nchooser *= n-k;
		nchooser /= k+1;
	}
	return nchooser;
}
inline void binary(UINT64 num)
{
	int rem;

	if (num <= 1)
	{
		std::cout<<num;
		return;
	}
	rem = num % 2;
	binary(num / 2);
	std::cout << rem;
}

inline UINT64 reorderbits(UINT64 bitstr1, UINT64 bitstr2, UINT64 chunk,int b)
{
	// b: total number of bits
	int tmp1, tmp2;
	UINT64 bitstr=0;
	for(int i=0;i<b;i++)
	{
		tmp1= bitstr1%2;
		tmp2 = bitstr2%2;
		if(chunk%2==0)
		{
			bitstr |= (UINT64)tmp1 << i;
			bitstr1 = bitstr1 >> 1;
		}

		if(chunk%2==1)
		{
			bitstr |= (UINT64)tmp2 << i;
			bitstr2= bitstr2 >> 1;
		}
		chunk = chunk>>1;

	}
	return bitstr;
}

inline int weight(UINT8 code){
    std::bitset<8> x;
    x = int (code);
    int wt=0;
    for (int i=0;i<8;i++){
        if (x[i]==1)
            wt+= (i%4)+1;
    }
    return wt;
}
inline void norm_chunks(UINT8* subnorms,UINT32 depth,UINT8* code,int B_over_8) {
	int numchunks = pow(2,depth);
	int Bchunks_over_8 = B_over_8/numchunks;
	int rem,div;
	UINT8 temp = 0;
	//if(depth==1)
	//	printf("as\n");
	//printf("here\n");   
    // numchunks=8  B_over_8 = 16 Bchunks_over_8=2
	if(numchunks<=B_over_8){
		for(int i =0;i<B_over_8;i++){
			rem = i%Bchunks_over_8;
			div = i/Bchunks_over_8;
            //std::cout<<weight((UINT32)code[i]) << " | "<<(UINT32)code[i]<<std::endl;
			temp += weight((UINT32)code[i]);
			if(rem == (Bchunks_over_8-1)){
				subnorms[div]=temp;
				temp=0;
			}
		}
	}
	else{
		UINT32 chunks_per_byte= numchunks/B_over_8;
		UINT8 mask;
		switch(chunks_per_byte) {
		case 2:
			mask = 0xF0; // mask = 11110000
			//printf("2 chunks per byte\n");
			break;
		case 4:
			mask = 0xC0; // mask = 11000000
			break;
		case 8:
			mask = 0x80; // mask = 10000000
			break;
		default:
			mask = 0;
			printf("Not Supported\n");
			return;
		}
		int index = 0;
		UINT8 masked;
		int bits_per_chunk = 8/chunks_per_byte;
		for(int i=0;i<B_over_8;i++){
			UINT8 smask = mask;
			for(UINT32 j=0;j<chunks_per_byte;j++)
			{
				masked = smask & code[i];
				subnorms[index] = weight((UINT32) masked);
				smask = smask >> bits_per_chunk;
				index++;
			}
		}


	}


}

inline int l1norm(UINT8* P,int codelb) {

	switch(codelb) {
	case 4: // 32 bit
		return popcnt(*(UINT32*)P);
	case 8: // 64 bit
		return popcntll(((UINT64*)P)[0]);
	case 16: // 128 bit
        std::cout<<popcntll(((UINT64*)P)[0]) \
				+ popcntll(((UINT64*)P)[1]) <<std::endl;
		return popcntll(((UINT64*)P)[0]) \
				+ popcntll(((UINT64*)P)[1]) ;
	case 32: // 256 bit
		return popcntll(((UINT64*)P)[0]) \
				+ popcntll(((UINT64*)P)[1])  \
				+ popcntll(((UINT64*)P)[2]) \
				+ popcntll(((UINT64*)P)[3]) ;
	case 64: // 512 bit
		return popcntll(((UINT64*)P)[0]) \
				+ popcntll(((UINT64*)P)[1])  \
				+ popcntll(((UINT64*)P)[2])  \
				+ popcntll(((UINT64*)P)[3])  \
				+ popcntll(((UINT64*)P)[4]) \
				+ popcntll(((UINT64*)P)[5]) \
				+ popcntll(((UINT64*)P)[6])  \
				+ popcntll(((UINT64*)P)[7]);
	default:
		int output = 0;
		for (int i=0; i<codelb; i++)
			output+= lookup[P[i]];
		return output;
	}

}

inline UINT64 pow2(int p) {
	UINT64 t = 1;
	return t<<p;
}




#endif /* BITOPS_H_ */
