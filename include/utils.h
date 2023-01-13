#ifndef UTILS_H_
#define UTILS_H_
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/time.h>
#include <cstring>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <random>
#include "../include/user_types.h"

void print(int);
float diff_ms(struct timespec *, struct timespec *);
template <typename T>
std::vector<T> to_vector(size_t, T *);
std::vector<float> to_vectorf32(int, float *);
std::vector<double> to_vectorf64(int, double *);
float_t *_read(std::string &, size_t, FILE *);
void do_fetch_data(raw_data_t *, point_t *, const size_t &);
#endif // UTILS_H_