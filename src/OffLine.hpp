#pragma once
// #ifndef __OFFLINE_H__
// #define __OFFLINE_H__

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
#include <smmintrin.h> 		// sse4.2
#include <immintrin.h>   	// avx
#include <pthread.h>

using namespace std;

#define unrolltime 2	//the times unrolling 
#define THREADNUM 4

//-------------The global variables for the seq and normal sse versions-------------------

// input params 
extern int 	state_num;	// number of states
extern int 	symbol_num;	// number of symbols
extern int 	enum_state;	// number of enumerative states (equal to symbol_num)
extern int 	spec_split1;	// number of species for spec version
extern int 	spec_split2;	// number of species for spec_avx version
extern int   	start_state;	// the start state
extern long 	len;  		// input array length
extern int 	KIND;		// the mode chosen

// runtime data 
extern int*	input __attribute__ ((aligned (32)));

// transition table for seq version
extern int*	T1 __attribute__ ((aligned (32)));  // dfa table 1d-state
extern int*	T2 __attribute__ ((aligned (32)));  // dfa table 1d-symbol
extern int**	T3 __attribute__ ((aligned (32)));  // dfa table 2d-state
extern int**	T4 __attribute__ ((aligned (32)));  // dfa table 2d-symbol

// parameter for spec
extern int*	predict_state;
extern int* spec_final;

// parameters for sse
extern uint8_t** 	T_dna __attribute__ ((aligned (8)));
extern uint8_t***	T_s;
extern int** 		U; // The unique state arrays  
extern uint8_t** 	U_bit;
extern int** 		L; // The look_up arrays 
extern int* 		NU; // The number of unique states for each column
extern int* 		NUV; //The correspoding 16-bit register used
extern int 		MAX_R;

// pthread version parameter
extern long charlen;
extern int* pthread_predict;
extern int* pthread_final;

//-------------Some instructions dealing with global variables for seq version------------
// load for T1 
void load_dfa_table1(char* dfafile1);

// load for T3
void load_dfa_table3(char* dfafile1);

// load for T4
void load_dfa_table4(char* dfafile1);

// load input to memory 
void load_input(char* inputfile) ;

// Constructing the dfa table with elements are 8-bits large 
void get_dfa();

//--------Instructions for completing the variables used in conv and randcoales versions-----------
	
// The offline work for range coalescing 
void OfflineWork();

//--------------SPEC version Offline------------------------

//-------------------------------------------
extern int* input_l; 
void inputlayout(int* a, int chunk);
void inputlayout_avxunroll(int* a, int chunk);


//--------------SPEC version Offline------------------------

void pthread_predict_op(int chunk);

void pthread_recompute(int chunk);

//#endif
