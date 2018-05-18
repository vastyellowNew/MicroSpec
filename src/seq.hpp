#pragma once
/* 
 * This files implements four kinds of seq version, including 1d-state.
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

#include "OffLine.hpp"

//1d-state
void seq();	
//1d-state for multi-threads
void* seq1(void *threadid);