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
#include <smmintrin.h> // sse4.2
#include <immintrin.h>   // avx

#include "OffLine.hpp"

// spec avx 1d-state version 
void* spec_avx1(void *threadid);
