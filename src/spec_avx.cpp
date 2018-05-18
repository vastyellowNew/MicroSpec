#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <cstdlib>
#include <sys/timeb.h>
#include <ctime>
#include <climits>
#include <iostream>
#include <fstream>
#include <vector>  
#include <algorithm>
#include <smmintrin.h> // sse4.2
#include <immintrin.h>   // avx

#include "spec_avx.hpp"

using namespace std;

// spec avx 1d-state version 
void* spec_avx1(void *threadid) 
{
	long tid =  (long)threadid;
	long i,j,k;

	__m256i scurrent_v;
	__m256i avxsymbol_v;
	__m256i avxmask_v;
	int* scurrent __attribute__ ((aligned (32)));
	scurrent= new int [spec_split2];
	int* current __attribute__ ((aligned (32)));
	current= new int [spec_split2];
	int* avxsymbol __attribute__ ((aligned (32)));
	avxsymbol=new int [spec_split2] ;
	int* addr __attribute__ ((aligned (32)));
	addr=new int [spec_split2];


	avxmask_v = _mm256_set1_epi32(0XFFFFFFFF);
	for(i=0;i<spec_split2;i++)
		scurrent[i]=pthread_predict[i+tid*spec_split2];

	scurrent_v=_mm256_maskload_epi32 ((int*)scurrent, avxmask_v);

	long bound=charlen/spec_split2;

	long cal[8];
	for(i=0;i<8;i++)
		cal[i]=i*bound;

	for(i = tid*charlen; i < bound+tid*charlen; i++)
	{
        avxsymbol_v=_mm256_set_epi32(input[i+cal[7]],input[i+cal[6]],input[i+cal[5]],input[i+cal[4]],input[i+cal[3]],input[i+cal[2]],input[i+cal[1]],input[i]);

        __m256i addr_v = _mm256_add_epi32(_mm256_mullo_epi32(scurrent_v, _mm256_set1_epi32(symbol_num)), avxsymbol_v);

		scurrent_v = _mm256_i32gather_epi32 (T1, addr_v, 4);

    }
    
    _mm256_store_si256((__m256i*)scurrent, scurrent_v);

    for(i = 0; i < spec_split2; i++)
		pthread_final[i+tid*spec_split2]=scurrent[i];

	printf("%ld is running on CPU %d\n", tid, sched_getcpu());
	pthread_exit((void*)threadid);
	
	delete []scurrent;
	delete []current;
	delete []avxsymbol;
	delete []addr;
}
