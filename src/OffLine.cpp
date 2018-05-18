
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

#include "OffLine.hpp"

using namespace std;

int 	state_num;	// number of states
int 	symbol_num;	// number of symbols
int 	enum_state;	// number of enumerative states (equal to symbol_num)
int 	spec_split1;	// number of species for spec version
int 	spec_split2;	// number of species for spec_avx version
int   	start_state;	// the start state
long 	len;  		// input array length
int 	KIND;		// the mode chosen

// runtime data 
int*	input __attribute__ ((aligned (32)));

// transition table for seq version
int*	T1 __attribute__ ((aligned (32)));  // dfa table 1d-state
int*	T2 __attribute__ ((aligned (32)));  // dfa table 1d-symbol
int**	T3 __attribute__ ((aligned (32)));  // dfa table 2d-state
int**	T4 __attribute__ ((aligned (32)));  // dfa table 2d-symbol

// parameter for spec
int*	predict_state;
int* spec_final;

// parameters for sse
uint8_t** 	T_dna __attribute__ ((aligned (8)));
uint8_t***	T_s;
int** 		U; // The unique state arrays  
uint8_t** 	U_bit;
int** 		L; // The look_up arrays 
int* 		NU; // The number of unique states for each column
int* 		NUV; //The correspoding 16-bit register used
int 		MAX_R=0;

// pthread version parameter
long charlen;
int* pthread_predict;
int* pthread_final;

//-------------Some instructions dealing with global variables for seq version------------
// load for T1 
void load_dfa_table1(char* dfafile1)
{
	int i,j;
	ifstream dfafile;
	dfafile.open(dfafile1);

	T1= new int [state_num*symbol_num];

	for(i=0;i<state_num;i++)
		for(j=0;j<symbol_num;j++)
		{
			int temp;
			dfafile>>temp;
			if(temp==(-1))
				temp=start_state;
			T1[i*symbol_num+j]=temp;
		}
	dfafile.close();
}

// load for T3
void load_dfa_table3(char* dfafile1)
{
	int i,j;
	ifstream dfafile;
	dfafile.open(dfafile1);

	T3 = new int* [state_num];
	for(i=0;i<state_num;i++)
		T3[i]=new int [symbol_num];

	for(i=0;i<state_num;i++)
		for(j=0;j<symbol_num;j++)
		{
			int temp;
			dfafile>>temp;
			if(temp==(-1))
				temp=start_state;
			T3[i][j]=temp;
		}								
	dfafile.close();
}

// load for T4
void load_dfa_table4(char* dfafile1)
{
	int i,j;
	ifstream dfafile;
	dfafile.open(dfafile1);

	T4 = new int* [symbol_num];
	for(i=0;i<symbol_num;i++)
		T4[i]=new int [state_num];

	for(i=0;i<state_num;i++)
		for(j=0;j<symbol_num;j++)
		{
			int temp;
			dfafile>>temp;
			if(temp==(-1))
				temp=start_state;
			T4[j][i]=temp;
		}
	dfafile.close();
}

// load input to memory 
void load_input(char* inputfile) 
{
	ifstream in;
	string str;
	long i=0;

	in.open(inputfile);
	len=0;
	while(getline(in,str))
	{
		len=str.size()+len;
	}
	in.close();

	in.open(inputfile);
	input = new int [len];	
	while(in)
	{
		char chara;
		in>>chara;

		//---------------For different benchmark-----------------------
		//DIV
		if(KIND==0)
		{
			if( chara == '0' )
				input[i] = 0;
			else if ( chara == '1' )
				input[i] = 1;
		}
		// SNORT
		else if(KIND==1)
		{
			input[i]=int(chara);
			if(input[i]>255||input[i]<0)
				i--;
		}
		//DNA
		else if(KIND==2)
		{
			if( chara == 'A' )
				input[i] = 0;
			else if ( chara == 'T' )
				input[i] = 1;
			else if ( chara == 'C' )
				input[i] = 2;
			else if ( chara == 'G' )
				input[i] = 3;
		}		
		// PROTN
		else if(KIND==3)
		{
			if( chara == 'A' )
				input[i] = 0;
			else if ( chara == 'C' )
				input[i] = 1;
			else if ( chara == 'D' )
				input[i] = 2;
			else if ( chara == 'E' )
				input[i] = 3;
			else if ( chara == 'F' )
				input[i] = 4;
			else if ( chara == 'G' )
				input[i] = 5;
			else if ( chara == 'H' )
				input[i] = 6;
			else if ( chara == 'I' )
				input[i] = 7;
			else if ( chara == 'K' )
				input[i] = 8;
			else if ( chara == 'L' )
				input[i] = 9;
			else if ( chara == 'M' )
				input[i] = 10;
			else if ( chara == 'N' )
				input[i] = 11;
			else if ( chara == 'P' )
				input[i] = 12;
			else if ( chara == 'Q' )
				input[i] = 13;
			else if ( chara == 'R' )
				input[i] = 14;
			else if ( chara == 'S' )
				input[i] = 15;
			else if ( chara == 'T' )
				input[i] = 16;
			else if ( chara == 'V' )
				input[i] = 17;
			else if ( chara == 'W' )
				input[i] = 18;
			else if ( chara == 'Y' )
				input[i] = 19;
		}
		// EVENODD
		else if (KIND==4)
		{
			if( chara == 'a' )
				input[i] = 0;
			else if ( chara == 'b' )
				input[i] = 1;
			else if ( chara == 'c' )
				input[i] = 2;
			else if ( chara == 'd' )
				input[i] = 3;
		}
		//---------------------------------------------------
		i++;
	}
	in.close();
}

// Constructing the dfa table with elements are 8-bits large 
void get_dfa()
{
	int i,j;
	T_dna=(uint8_t **)malloc(sizeof(uint8_t *) * symbol_num);
    for(i = 0; i < symbol_num; i++)
        *(T_dna + i) = (uint8_t *)malloc(sizeof(uint8_t) * state_num); 
	for(i=0;i<symbol_num;i++)
        	for(j=0;j < state_num;j++)
               		T_dna[i][j]=T4[i][j];
}

//--------Instructions for completing the variables used in conv and randcoales versions-----------
	
// The offline work for range coalescing 
void OfflineWork()
{	
	cout<<"Offline work."<<endl; 
	int i,j,k,l;
	int j1,j2,j3;
	vector <int> vec;

	// Getting the Unique arrays and look_up arrays
	U=new int* [symbol_num];
	U_bit= (uint8_t**)malloc(sizeof(uint8_t*) * symbol_num);
	L=new int* [symbol_num];
	NU=new int [symbol_num];
	NUV=new int [symbol_num];

	for(i=0;i<symbol_num;i++)
	{
		for(j1=0;j1<state_num;j1++)	
			if(find(vec.begin(),vec.end(),T4[i][j1])==vec.end())
				vec.push_back(T4[i][j1]);
		int sizeofu;
		sizeofu=vec.size();
		U[i]=new int [sizeofu];
		*(U_bit+i)=(uint8_t*)malloc(sizeof(uint8_t) * sizeofu);
		L[i]=new int [state_num];
		NU[i]=sizeofu;
		
		if(MAX_R<sizeofu)	
			MAX_R=sizeofu;

		if(NU[i]%16>0)
			NUV[i]=NU[i]/16+1;
		else NUV[i]=NU[i]/16;
		
		for(j2=0;j2<sizeofu;j2++)
		{
			U[i][j2]=vec[j2];
			cout<<j2<<endl;
		}

		for(j3=0;j3<state_num;j3++)
			for(k=0;k<sizeofu;k++)
				if(T4[i][j3]==U[i][k])
					L[i][j3]=U[i][k];
		vec.erase(vec.begin(),vec.end());
	}		

	// Constructing the small transition table 
	T_s=(uint8_t ***)malloc(sizeof(uint8_t **) * symbol_num);
    	for(i = 0; i < symbol_num; i++)
         	T_s[i] = (uint8_t **)malloc(sizeof(uint8_t*) * symbol_num); 
	for(i = 0; i < symbol_num; i++)
		for(j = 0; j < symbol_num; j++)
			T_s[i][j]= (uint8_t *)malloc(sizeof(uint8_t) * NU[i]);

    	for(i = 0; i < symbol_num; i++)
		for(j = 0; j < symbol_num; j++)
			for(k = 0; k < NU[i]; k++)
				for(l = 0; l < NU[j]; l++)
					if(T4[j][U[i][k]]==U[j][l])
						T_s[i][j][k]=(uint8_t)l;
	cout<<"Max range is "<<MAX_R<<endl;
}


//-----------------Trans----------------

int* input_l; 
void inputlayout(int* a, int chunk)
{
	input_l=new int [len];
	long tcharlen=len/THREADNUM;
	long bound=tcharlen/chunk;

	for (int tc=0;tc<THREADNUM;tc++)
	{
		for(long i=0;i<bound;i++)
			for(int j=0;j<chunk;j++)
				input_l[i*chunk+j+tc*tcharlen]=a[j*bound+i+tc*tcharlen];
	}
	//delete []a;
	//return b;
}

void inputlayout_avxunroll(int* a, int chunk)
{
	input_l=new int [len];
	long tcharlen=len/THREADNUM;
	long bound=tcharlen/chunk;
	
	for (int tc=0;tc<THREADNUM;tc++)
	{
		for(long i=0;i<bound;i++)
			for(int j=0;j<chunk;j++)
	                {
				if(j%2==0)
					input_l[i*chunk+j/2+tc*tcharlen]=a[j*bound+i+tc*tcharlen];
		                else
		                        input_l[i*chunk+chunk/2+(j-1)/2+tc*tcharlen]=a[j*bound+i+tc*tcharlen];
		         }
	}
	//delete []a;
	//return b;
}


//--------------SPEC version Offline------------------------


//--------------SPEC version Offline------------------------

void pthread_predict_op(int chunk)
{
	int i,j;
	pthread_predict=new int [chunk];
	long pcharlen=len/chunk;

	pthread_predict[0]=start_state;
	for (j=1;j<chunk;j++)
	{
		int look=len*j/chunk;
		int ite=start_state;
		for(i=0;i<1000;i++)
		{
			int symbol=input[look-1000+i];
			ite=T3[ite][symbol];
		}
		pthread_predict[j]=ite;
	}
	for(j=0;j<chunk;j++)
		cout<<"The start state of chunk "<<j<<" is "<<pthread_predict[j]<<endl;
}

void pthread_recompute(int chunk)
{
	int c1;
	for (c1=0;c1<chunk;c1++)
		cout<<"The final state of "<<c1<<" thread is "<<pthread_final[c1]<<endl;

	long i;
	long pcharlen;
	pcharlen=len/chunk;
	for(c1=0;c1<chunk-1;c1++)
	{
		if(pthread_final[c1] != pthread_predict[c1+1])
		{
			i=pcharlen*(c1+1);
			int temp1=pthread_final[c1];
			int temp2=pthread_predict[c1+1];
			do
			{
				temp1=T3[temp1][input[i]];
				temp2=T3[temp2][input[i]];
				i++;
				if(temp1==temp2)
					break;
			}while(i<pcharlen*(c1+2));
			cout<<"The reprocessing length in chunk "<<c1+1<<" is "<<(i-pcharlen*(c1+1))<<endl;
		}
	}
}

//#endif
