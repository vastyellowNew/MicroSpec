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

#include "spec_avxunroll.hpp"

using namespace std;

// Speculation avx vectorization version wiht unrolling
void* spec_avx_unroll1(void *threadid) 
{
	long tid =  (long)threadid;
	long i,j ;
	int k;

	__m256i current_v[unrolltime];
	__m256i symbol_v[unrolltime];
    __m256i addr_v[unrolltime];
	__m256i mask;

	int* current __attribute__ ((aligned (32)));
	current= new int [spec_split2];
	int* scurrent __attribute__ ((aligned (32)));
	scurrent= new int [spec_split2];
	int* symbol __attribute__ ((aligned (32)));
	symbol=new int [spec_split2] ;
	int* addr __attribute__ ((aligned (32)));
	addr=new int [spec_split2];

	mask = _mm256_set1_epi32(0XFFFFFFFF);
	current[0]=start_state;

	
	for(k=0;k<unrolltime;k++)
	{
		for(i=0;i<spec_split2;i++)
			current[i]=pthread_predict[i*unrolltime+tid*spec_split2*unrolltime+k];
		current_v[k]=_mm256_maskload_epi32 ((int*)current, mask);
	}

	long bound = charlen/(unrolltime*spec_split2);
	int cal[unrolltime][8];
	for(i=0;i<unrolltime;i++)
		for(j=0;j<8;j++)
			cal[i][j]=i*bound+j*(charlen/spec_split2);

	for(i = tid*charlen; i < bound+tid*charlen; i++)
	{
        	for(k=0;k<unrolltime;k++)
        	{
        		symbol_v[k]=_mm256_set_epi32(input[i+cal[k][7]],input[i+cal[k][6]],input[i+cal[k][5]],input[i+cal[k][4]],input[i+cal[k][3]],input[i+cal[k][2]],input[i+cal[k][1]],input[i+cal[k][0]]);

            	addr_v[k] = _mm256_add_epi32(_mm256_mullo_epi32(current_v[k], _mm256_set1_epi32(symbol_num)), symbol_v[k]);
            	current_v[k] = _mm256_i32gather_epi32 (T1, addr_v[k], 4);
         	}
    }

    for(k=0;k<unrolltime;k++)
    {
        _mm256_store_si256((__m256i*)current, current_v[k]);
        for(i = 0; i < spec_split2; i++)
		{        		
			pthread_final[k+i*unrolltime+tid*spec_split2*unrolltime]=current[i];
		}
    }
    
    printf("%ld is running on CPU %d\n", tid, sched_getcpu());
	pthread_exit((void*)threadid);

	delete []current;
	delete []scurrent;
	delete []symbol;
	delete []addr;
}

//----------------------level4 unroll and avx-----------------------------------------
// Speculation avx vectorization version wiht unrolling
void* spec_unroll_avx1(void *threadid) 
{
	long tid =  (long)threadid;
	long i,j ;
	int k;

	__m256i current_v[unrolltime];
	__m256i symbol_v[unrolltime];
    __m256i addr_v[unrolltime];
	__m256i mask;

	int* current __attribute__ ((aligned (32)));
	current= new int [spec_split2];
	int* scurrent __attribute__ ((aligned (32)));
	scurrent= new int [spec_split2];
	int* symbol __attribute__ ((aligned (32)));
	symbol=new int [spec_split2] ;
	int* addr __attribute__ ((aligned (32)));
	addr=new int [spec_split2];

	mask = _mm256_set1_epi32(0XFFFFFFFF);
	
	for(k=0;k<unrolltime;k++)
	{
		for(i=0;i<spec_split2;i++)
			current[i]=pthread_predict[i+tid*spec_split2*unrolltime+k*spec_split2];
		current_v[k]=_mm256_maskload_epi32 ((int*)current, mask);
	}
		

	long bound = charlen/(unrolltime*spec_split2);
	long wall=charlen/unrolltime;
	int cal[unrolltime][8];
	for(i=0;i<unrolltime;i++)
		for(j=0;j<8;j++)
			cal[i][j]=i*wall+j*bound;
	
	for(i = tid*charlen; i < bound+(tid)*charlen; i++)
	{
        	for(k=0;k<unrolltime;k++)
        	{
        		symbol_v[k]=_mm256_set_epi32(input[i+cal[k][7]],input[i+cal[k][6]],input[i+cal[k][5]],input[i+cal[k][4]],input[i+cal[k][3]],input[i+cal[k][2]],input[i+cal[k][1]],input[i+cal[k][0]]);

            	addr_v[k] = _mm256_add_epi32(_mm256_mullo_epi32(current_v[k], _mm256_set1_epi32(symbol_num)), symbol_v[k]);
            	current_v[k] = _mm256_i32gather_epi32 (T1, addr_v[k], 4);
         	}
    }
    
    for(k=0;k<unrolltime;k++)
    {
        _mm256_store_si256((__m256i*)current, current_v[k]);
        for(i = 0; i < spec_split2; i++)
		{
			pthread_final[i+k*spec_split2+tid*spec_split2*unrolltime]=current[i];
		}
    }

    printf("%ld is running on CPU %d\n", tid, sched_getcpu());
	pthread_exit((void*)threadid);

	delete []current;
	delete []scurrent;
	delete []symbol;
	delete []addr;
}
