#ifndef INTERPOLATION_H__
#define INTERPOLATION_H__
#include "user_types.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

void parabolic(FFT_data *, int);
void gaussian(FFT_data *, int);
float Neville(float *, float *, float, int);
#endif
