#ifndef USER_TYPES_H__
#define USER_TYPES_H__
#include <math.h>
#include <string.h>
#pragma once
#define buffer_size 1000

typedef struct raw_data_t
{
    std::string path_t;
    std::string path_y;
    float_t *ys = NULL, *ts = NULL, avg = 0;
    size_t file_length;
} raw_data_t;

typedef struct freq_t
{
    size_t zc_counts, zc_indices[buffer_size] = {0};
    float_t
        *freq,
        mean_freq = 0,
        std_freq = 0,
        prec_in_t = 0,
        freqs[buffer_size] = {0},
        zc_ts[buffer_size] = {0},
        zc_ys[buffer_size] = {0};
} freq_t;

typedef struct point_t
{
    float_t *ys = NULL, *ts = NULL;
    bool TRANSFER_CHECK = false;
    size_t length = 0;
} point_t;
#endif