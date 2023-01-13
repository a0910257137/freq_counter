#ifndef WAVE_GENERATOR_H__
#define WAVE_GENERATOR_H__
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <iostream>

typedef struct wave_t
{
    float_t *ys = NULL, *ts = NULL, average = 0;
    size_t arr_size;
} wave_t;

double_t f(double_t);
wave_t *gen_waves(float_t);
#endif