#ifndef FFT_H__
#define FFT_H__
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "utils.h"
#include "user_types.h"
#include <sys/stat.h>
#include <sys/types.h>
#ifndef PI
#define PI 3.14159265358979323846264338327950288
#endif

void _get_common_table(float *, unsigned int &, unsigned int);
void improve_fft(float *, float *, float *, float *, float *, unsigned int &, unsigned int);
void fft(complex *, int, complex *);
void _peaks(FFT_data *, int, int);
void __init_fft_params(FFT_data *, const unsigned int &, const unsigned int &);
void __reset_params(FFT_data *, const unsigned int &);
void __free_params(FFT_data *);
#endif