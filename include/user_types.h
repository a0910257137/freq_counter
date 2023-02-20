#ifndef USER_TYPES_H__
#define USER_TYPES_H__
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#pragma once
#define buffer_size 1000

typedef float real;
typedef struct
{
    real Re;
    real Im;
} complex;

typedef struct
{
    complex *fft_vals;
    float *Re, *Im, *temp_Re, *temp_Im, *pre_build_table;
    float *k_peaks, *k_frequencies, *k_maginitudes;
    float *magnitudes, *freq_steps, peak_magnitude;
    float T, estimated_freq;
    unsigned int peak, top_k, table_idx;
} FFT_data;

typedef struct zc_t
{
    int zc_counts, *zc_indices;
    unsigned int K, inter_N;
    float
        *freqs,
        *zc_ts,
        *zc_ys,
        mean_freq,
        std_freq,
        time_step;
} zc_t;

#endif