#ifndef ZERO_CROSSING_H__
#define ZERO_CROSSING_H__
#include <string.h>
#include <math.h>
#include "user_types.h"
#include "interpolation.h"
#include "utils.h"
#define buffer_size 1000
int sgn(float);
void estimate(zc_t *, float *, const unsigned int &);
void __init_zc_params(zc_t *, const unsigned int &);
void __reset_zc_params(zc_t *src_data, const unsigned int &Fs);
#endif