/*
	* This file can be applied to the general benchmarks. 

	* The inputs you have to offer are 
	* Dfa file name, Input file name, 
	* number of states, number of symbols, number of species for speculation algorithm, 
	* start state for seq version, the kind of benchmark, and the mode.
	* This version is implemented for using 4 threads (It can be modified for more threads).
	* edited by J.Q. All copyright reserved.
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
#include <math.h>
#include <smmintrin.h> 		// sse4.2
#include <immintrin.h>   	// avx
#include <pthread.h>

#include "OffLine.hpp"
#include "seq.hpp"
#include "spec_nonSIMD.hpp"
#include "spec_avx.hpp"
#include "spec_avxunroll.hpp"
#include "sse.hpp"

using namespace std;	

#define THREADNUM 4	

int main(int argv, char* argc[])
{
	char* dfafile;
	char* inputfile;
	ofstream outfile1,outfile2;

	//Initialization
	dfafile=argc[1]; 			// The dfa transtion file name
	inputfile=argc[2]; 			// The input file name
	state_num=atoi(argc[3]);	// The number of states
	symbol_num=atoi(argc[4]);	// The number of symbols
	spec_split1=atoi(argc[5]);	// The number of species (for normal speculation)
	spec_split2=atoi(argc[6]); 	// The number of species (for avx speculation)
	start_state=atoi(argc[7]);	// The start state
	KIND=atoi(argc[8]);			// The kind of benchmark
	int mode=atoi(argc[9]);		// The mode
	enum_state=state_num;		// The number of enumerative states in sse version
	
	//Load operation
    load_input(inputfile);
	load_dfa_table3(dfafile);	// Loading 2D state-major table
	load_dfa_table4(dfafile);	// Loading 2D symbol-major table

	// In this experiments, ensure that the input can be divided by spec_split	
	//len=100000000;
	len=(len/(spec_split1*spec_split2))*spec_split1*spec_split2; 
    cout<<endl<<"The input length is "<<len<<endl<<endl;
    long cmp1, cmp1r;
    struct timeb startTime1, endTime1;
    struct timeb startTime1r, endTime1r;

    // pthread version variables
	int rc;
   	int temptb;
   	long pt;
   	pthread_t* threads;
	cpu_set_t* cpu;

	threads=(pthread_t*)malloc(sizeof(pthread_t)*THREADNUM);	// thread set;
	cpu=(cpu_set_t*)malloc(sizeof(cpu_set_t)*THREADNUM);		// thread binding variables
	charlen = len/THREADNUM;

	// thread binding
	for(pt=0; pt<THREADNUM; pt++)
	{
    		CPU_ZERO(&cpu[pt]);
    		CPU_SET(pt, &cpu[pt]);
	}

//----------seq version------------------------------------------
	if (mode==0)
	{		
		load_dfa_table1(dfafile);
		// Time Spent
		ftime(&startTime1);
		seq();		
		ftime(&endTime1);

		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cout<<cmp1<<endl;
		outfile1.open("R_seq.txt", ios::app);
		outfile1 << cmp1 << endl;
		outfile1.close();

		delete []T1;
		cout<<endl;
	}

//----------Multi-threads seq version------------------------------------------
	if (mode==1)
	{		
		load_dfa_table1(dfafile);
		pthread_predict_op(THREADNUM);
		pthread_final=new int [THREADNUM];

		// Time Spent
		ftime(&startTime1);
		for(pt=0;pt<THREADNUM;pt++)
		{
     			printf("In main: creating thread %ld\n", pt);
     			rc = pthread_create(&threads[pt], NULL, seq1, (void *)pt);
     			if (rc)
				{
       				printf("ERROR; return code from pthread_create() is %d\n", rc);
       				exit(-1);
       			}
				// thread binding
    			temptb = pthread_setaffinity_np(threads[pt], sizeof(cpu_set_t), &cpu[pt]);
     	}
		for(pt=0; pt<THREADNUM; pt++)
        		pthread_join(threads[pt], NULL);	
		ftime(&endTime1);

		ftime(&startTime1r);
		pthread_recompute(THREADNUM);
		ftime(&endTime1r);

		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cmp1r=(endTime1r.time-startTime1r.time)*1000+(endTime1r.millitm-startTime1r.millitm);
		cout<<cmp1+cmp1r<<endl;
		cout<<"Re time "<<cmp1r<<endl;
		outfile1.open("R_seq_MultiT.txt", ios::app);
		outfile1 << cmp1+cmp1r << endl;
		outfile1.close();

		delete []T1;
		delete []pthread_predict;
		delete []pthread_final;
		cout<<endl;
	}

//---------------------spec nonSIMD version-------------------------------
	if(mode==2)
	{
		load_dfa_table1(dfafile);
		pthread_predict_op(THREADNUM*spec_split1);
		pthread_final=new int [THREADNUM*spec_split1];
		
		ftime(&startTime1);
		
		for(pt=0;pt<THREADNUM;pt++)
		{
     			printf("In main: creating thread %ld\n", pt);
     			rc = pthread_create(&threads[pt], NULL, spec1, (void *)pt);
     			if (rc)
				{
       				printf("ERROR; return code from pthread_create() is %d\n", rc);
       				exit(-1);
       			}
				// thread binding
    			temptb = pthread_setaffinity_np(threads[pt], sizeof(cpu_set_t), &cpu[pt]);
     	}
		for(pt=0; pt<THREADNUM; pt++)
        	pthread_join(threads[pt], NULL);

		ftime(&endTime1);
		
		ftime(&startTime1r);
		pthread_recompute(THREADNUM*spec_split1);
		ftime(&endTime1r);

		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cmp1r=(endTime1r.time-startTime1r.time)*1000+(endTime1r.millitm-startTime1r.millitm);
		cout<<cmp1+cmp1r<<endl;
		cout<<"Re time "<<cmp1r<<endl;
		
		outfile1.open("R_unroll_multi.txt", ios::app);
		outfile1 << cmp1+cmp1r << endl;
		outfile1.close();

		delete []T1;
		delete []pthread_predict;
		delete []pthread_final;
		cout<<endl;
	}

//-------------------------spec avx-------------------------------
	if(mode==3)
	{
		load_dfa_table1(dfafile);
		pthread_predict_op(THREADNUM*spec_split2);
		pthread_final=new int [THREADNUM*spec_split2];

		ftime(&startTime1);

		for(pt=0;pt<THREADNUM;pt++)
		{
     			printf("In main: creating thread %ld\n", pt);
     			rc = pthread_create(&threads[pt], NULL, spec_avx1, (void *)pt);
     			if (rc)
				{
       				printf("ERROR; return code from pthread_create() is %d\n", rc);
       				exit(-1);
       			}
				// thread binding
    			temptb = pthread_setaffinity_np(threads[pt], sizeof(cpu_set_t), &cpu[pt]);
     	}
		for(pt=0; pt<THREADNUM; pt++)
        	pthread_join(threads[pt], NULL);

		ftime(&endTime1);
		
		ftime(&startTime1r);
		pthread_recompute(THREADNUM*spec_split2);
		ftime(&endTime1r);
		
		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cmp1r=(endTime1r.time-startTime1r.time)*1000+(endTime1r.millitm-startTime1r.millitm);
		cout<<cmp1+cmp1r<<endl;
		cout<<"Re time "<<cmp1r<<endl;
		cout<<endl;
		outfile1.open("R_avx_multi.txt", ios::app);
		outfile1 << cmp1+cmp1r << endl;
		outfile1.close();

		delete []T1;
		delete []pthread_predict;
		delete []pthread_final;
	}

//-----------------------sse version-------------------------------
	// sse normal version
	if(mode==5)
	{
		get_dfa();

		ftime(&startTime1);

		for(pt=0;pt<THREADNUM;pt++)
		{
     			printf("In main: creating thread %ld\n", pt);
     			rc = pthread_create(&threads[pt], NULL, sse_dna, (void *)pt);
     			if (rc)
				{
       				printf("ERROR; return code from pthread_create() is %d\n", rc);
       				exit(-1);
       			}
				// thread binding
    			temptb = pthread_setaffinity_np(threads[pt], sizeof(cpu_set_t), &cpu[pt]);
     	}
		for(pt=0; pt<THREADNUM; pt++)
        	pthread_join(threads[pt], NULL);

		ftime(&endTime1);
		
		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cout<<cmp1<<endl;
		cout<<endl;
		outfile1.open("R_sse_normal.txt", ios::app);
		outfile1 << cmp1 << endl;
		outfile1.close();
 	}

	// sse with range coalescing version 
	if(mode==6)
	{
		get_dfa();
		OfflineWork();

		ftime(&startTime1);

		for(pt=0;pt<THREADNUM;pt++)
		{
     			printf("In main: creating thread %ld\n", pt);
     			rc = pthread_create(&threads[pt], NULL, sse_dna_rangecoalescing, (void *)pt);
     			if (rc)
				{
       				printf("ERROR; return code from pthread_create() is %d\n", rc);
       				exit(-1);
       			}
				// thread binding
    			temptb = pthread_setaffinity_np(threads[pt], sizeof(cpu_set_t), &cpu[pt]);
     	}
		for(pt=0; pt<THREADNUM; pt++)
        	pthread_join(threads[pt], NULL);
		
		ftime(&endTime1);
		
		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cout<<cmp1<<endl;
		
		outfile1.open("R_sse_range.txt", ios::app);
		outfile1 << cmp1 << endl;
		outfile1.close();

		cout<<endl;
	}
	
	// sse with convergence version  
	if(mode==7) 
	{
		get_dfa();

		ftime(&startTime1);

		for(pt=0;pt<THREADNUM;pt++)
		{
     			printf("In main: creating thread %ld\n", pt);
     			rc = pthread_create(&threads[pt], NULL, sse_dna_convergence, (void *)pt);
     			if (rc)
				{
       				printf("ERROR; return code from pthread_create() is %d\n", rc);
       				exit(-1);
       			}
				// thread binding
    			temptb = pthread_setaffinity_np(threads[pt], sizeof(cpu_set_t), &cpu[pt]);
     	}
		for(pt=0; pt<THREADNUM; pt++)
        	pthread_join(threads[pt], NULL);

		ftime(&endTime1);
		
		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cout<<cmp1<<endl;
		cout<<endl;
		
		outfile1.open("R_sse_converg.txt", ios::app);
		outfile1 << cmp1 << endl;
		outfile1.close();
	} 

//-----------------------avx+unroll version-------------------------------

	if(mode==8)
	{
		load_dfa_table1(dfafile);
		pthread_predict_op(THREADNUM*spec_split1*unrolltime);
		pthread_final=new int [THREADNUM*spec_split2*unrolltime];

		ftime(&startTime1);

		for(pt=0;pt<THREADNUM;pt++)
		{
     			printf("In main: creating thread %ld\n", pt);
     			rc = pthread_create(&threads[pt], NULL, spec_avx_unroll1, (void *)pt);
     			if (rc)
				{
       				printf("ERROR; return code from pthread_create() is %d\n", rc);
       				exit(-1);
       			}
				// thread binding
    			temptb = pthread_setaffinity_np(threads[pt], sizeof(cpu_set_t), &cpu[pt]);
     	}
		for(pt=0; pt<THREADNUM; pt++)
        	pthread_join(threads[pt], NULL);
		
		ftime(&endTime1);
		
		ftime(&startTime1r);
		pthread_recompute(THREADNUM*spec_split2*unrolltime);
		ftime(&endTime1r);

		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cmp1r=(endTime1r.time-startTime1r.time)*1000+(endTime1r.millitm-startTime1r.millitm);
		cout<<cmp1+cmp1r<<endl;
		cout<<"Re time "<<cmp1r<<endl;
		cout<<endl;

		outfile1.open("R_avxunroll_multi.txt", ios::app);
		outfile1 << cmp1+cmp1r << endl;
		outfile1.close();

		delete []T1;
		delete []pthread_predict;
		delete []pthread_final;
	}

	// Unrolling and AVX speculation version 1d-state
	if(mode==10)
	{

		load_dfa_table1(dfafile);
		pthread_predict_op(THREADNUM*spec_split2*unrolltime);
		pthread_final=new int [THREADNUM*spec_split2*unrolltime];

		ftime(&startTime1);
				
		for(pt=0;pt<THREADNUM;pt++)
		{
     			printf("In main: creating thread %ld\n", pt);
     			rc = pthread_create(&threads[pt], NULL, spec_unroll_avx1, (void *)pt);
     			if (rc)
				{
       				printf("ERROR; return code from pthread_create() is %d\n", rc);
       				exit(-1);
       			}
				// thread binding
    			temptb = pthread_setaffinity_np(threads[pt], sizeof(cpu_set_t), &cpu[pt]);
     	}
		for(pt=0; pt<THREADNUM; pt++)
        	pthread_join(threads[pt], NULL);
		
		ftime(&endTime1);
		
		ftime(&startTime1r);
		pthread_recompute(THREADNUM*spec_split2*unrolltime);
		ftime(&endTime1r);

		cmp1=(endTime1.time-startTime1.time)*1000+(endTime1.millitm-startTime1.millitm);
		cmp1r=(endTime1r.time-startTime1r.time)*1000+(endTime1r.millitm-startTime1r.millitm);
		cout<<cmp1+cmp1r<<endl;
		cout<<"Re time "<<cmp1r<<endl;
		cout<<endl;
		outfile1.open("R_unrollavx_multi.txt", ios::app);
		outfile1 << cmp1+cmp1r << endl;
		outfile1.close();

		delete []T1;
		delete []pthread_predict;
		delete []pthread_final;	
	}

	return 0;
}
