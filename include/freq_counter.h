#ifndef FREQ_COUNTER_H__
#define FREQ_COUNTER_H__
#include <string.h>
#include "../include/utils.h"
#include "matplotlibcpp.h"
#include <vector>
#include "../include/user_types.h"
template <typename T>
int sgn(T);
freq_t *estimate(point_t *, size_t, size_t);
#endif