#ifndef CMPLX_FUNCS_H__
#define CMPLX_FUNCS_H__
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "user_types.h"

void cmplx_peak_mag(float *, float *, float *, float &, unsigned int &, const unsigned int);
void find_k_peaks(FFT_data *, const unsigned int);
#endif