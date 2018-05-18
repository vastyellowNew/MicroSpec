/* 
 * This files implements four kinds of seq version, including 1d-state, 1d-symbol, 2d-state, 2d-symbol and 2d-symbol-shuffle seq.
 * The output operation is store operation, which means store the current state for each iteration. But here we do not concenned about it.
 * This file is edited by J.Q. All copyright reserved.
*/

#include <stdio.h>
#include <stdint.h>
#include <string>
#include <cstdlib>
#include <sys/timeb.h>
#include <ctime>
#include <climits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>  
#include <algorithm>
#include <smmintrin.h> // sse4.2
#include <immintrin.h>   // avx
#include <math.h>

#include "seq.hpp"

using namespace std;	

//1d-state
void seq()
{
    long i;
	int current,symbol;
	current=start_state;

	for(i=0;i<len;i++)
	{
		symbol=input[i];
		current=T1[current*symbol_num+symbol];
	}
	cout<<"Seq1 Final State is "<<current<<endl;
}	

//-----------------Different versions--------------------------
//1d-state for multi-threads
void* seq1(void *threadid)
{
	long tid =  (long)threadid;
    long i;
	int current,symbol;
	current=pthread_predict[tid];

	long pbound=(tid+1)*charlen;
	long starti=tid*charlen;
	for(i=starti;i<pbound;i++)
	{
		symbol=input[i];
		current=T1[current*symbol_num+symbol];
	}

	pthread_final[tid]=current;
    printf("%ld is running on CPU %d\n", tid, sched_getcpu());
	pthread_exit((void*)threadid);
}
