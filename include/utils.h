#ifndef UTILS_H_
#define UTILS_H_
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "user_types.h"
#include <vector>
extern void print(int);
extern float diff_ms(struct timespec *, struct timespec *);
template <typename T>
std::vector<T> to_vector(int, T *);
std::vector<float> to_vectorf32(int, float *);
std::vector<double> to_vectorf64(int, double *);
float *_read(const char *, int, FILE *);
void do_fetch_data(float *, float *, unsigned const int);
void min_max_norm(float *, int);
#endif // UTILS_H_