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

#include "spec_nonSIMD.hpp"

using namespace std;

// Normal speculation version
void* spec1(void *threadid)
{
	long tid =  (long)threadid;
	int* scurrent;
	int* avxsymbol;
	long i,j,gap;

	scurrent = new int [spec_split1];
	avxsymbol = new int [spec_split1];
	
	for(j=0;j<spec_split1;j++)
		scurrent[j]=pthread_predict[j+tid*spec_split1];
	
	long bound=charlen/spec_split1;
	for(i=tid*charlen;i<bound+tid*charlen;i++)
	{
		gap = 0;
		scurrent[0]=T1[scurrent[0]*symbol_num+input[i+gap]];
		gap = gap + bound;
		scurrent[1]=T1[scurrent[1]*symbol_num+input[i+gap]];
		gap = gap + bound;
		scurrent[2]=T1[scurrent[2]*symbol_num+input[i+gap]];
		gap = gap + bound;
		scurrent[3]=T1[scurrent[3]*symbol_num+input[i+gap]];
		gap = gap + bound;
		scurrent[4]=T1[scurrent[4]*symbol_num+input[i+gap]];
		gap = gap + bound;
		scurrent[5]=T1[scurrent[5]*symbol_num+input[i+gap]];
		gap = gap + bound;
		scurrent[6]=T1[scurrent[6]*symbol_num+input[i+gap]];
		gap = gap + bound;
		scurrent[7]=T1[scurrent[7]*symbol_num+input[i+gap]];

	}
	for(i=0;i<spec_split1;i++)
		pthread_final[tid*spec_split1+i]=scurrent[i];

	printf("%ld is running on CPU %d\n", tid, sched_getcpu());
	pthread_exit((void*)threadid);

	delete []scurrent;
	delete []avxsymbol;
}
