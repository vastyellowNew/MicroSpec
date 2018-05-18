#pragma once

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
#include <math.h>
#include <smmintrin.h> // sse4.2
#include <immintrin.h>   // avx

#include "OffLine.hpp"

void* sse_dna(void *threadid);
void* sse_dna_rangecoalescing(void* threadid);
void* sse_dna_convergence(void* threadid);
