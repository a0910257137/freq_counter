#include "matplotlibcpp.h"
#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include "../include/user_types.h"
#include "../include/freq_counter.h"

using namespace std;
using std::vector;
namespace plt = matplotlibcpp;
timespec ivk_s, ivk_e;
#define CHECK(x)                                                 \
    if (!(x))                                                    \
    {                                                            \
        fprintf(stderr, "Error at %s:%d\n", __FILE__, __LINE__); \
        exit(1);                                                 \
    }

int main(int argc, char *argv[])
{
    // assume coming data type is float 4bytes
    // **********************************
    //  frame = 32 [unit:MBytes]
    //  N = 1, 2, 3, 4, 5
    //**********************************
    raw_data_t *raw_data = (raw_data_t *)malloc(1 * sizeof(raw_data_t));
    FILE *ptr;
    raw_data->path_t = "/aidata/anders/objects/freqs/1M_250M_ts_rn.bin";
    raw_data->path_y = "/aidata/anders/objects/freqs/1M_250M_ys_rn.bin";
    raw_data->file_length = 250 * 1e6;
    float_t *ts = _read(raw_data->path_t, raw_data->file_length, ptr);
    float_t *ys = _read(raw_data->path_y, raw_data->file_length, ptr);

    raw_data->ts = ts;
    raw_data->ys = ys;

    const size_t N = 1, unit = 32 * 1e6;
    point_t *points = (point_t *)malloc(1 * sizeof(point_t));
    points->length = N * unit;
    points->ts = (float_t *)malloc(points->length * sizeof(float_t));
    points->ys = (float_t *)malloc(points->length * sizeof(float_t));
    int i = 0;
    while (1)
    {
        do_fetch_data(raw_data, points, points->length);
        //********************use the noise data to estimate ********************
        // vector<float>
        //     t_vec = to_vectorf32(points->length, points->ts);
        // vector<float>
        //     y_vec = to_vectorf32(points->length, points->ys);
        // plt::ylabel("unit: V"); // Plot the label
        // plt::xlabel("unit: s"); // Plot the label
        // plt::xlim(0.00, 2e-6);
        // // Set y-axis interval
        // plt::plot(t_vec, y_vec);
        // plt::legend();
        // plt::xlabel("Time");
        // plt::ylabel("Voltage");
        // plt::grid(true);
        // plt::save("/home2/anders/proj_c/CCounter/basic.png");
        // exit(1);
        clock_gettime(CLOCK_REALTIME, &ivk_s);
        freq_t *freq_counter = estimate(points, 5, 10);
        clock_gettime(CLOCK_REALTIME, &ivk_e);
        float diff = diff_ms(&ivk_e, &ivk_s);
        std::cout << "cost time: " << diff << (" ms") << std::endl;
        std::cout
            << "estimate frequency: " << freq_counter->mean_freq << " Hz" << std::endl;
        free(freq_counter);
    }
    free(raw_data);
}