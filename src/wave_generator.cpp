#include "../include/wave_generator.h"
double f(double x)
{
    if (x >= 0 && x < M_PI)
    {
        return x;
    }
    else if (x >= M_PI && x < 2 * M_PI)
    {
        return 2 * M_PI - x;
    }
    else if (x >= 2 * M_PI)
    {
        return f(x - 2 * M_PI);
    }
    else if (x < 0)
    {
        return f(x + 2 * M_PI);
    }
}
wave_t *gen_waves(float_t freq)
{

    // sampling rates = 1 M/sec
    // assuming only 0.01 seconds
    wave_t *outputs = (wave_t *)malloc(1 * sizeof(wave_t));
    float_t measured_sec = 0.1, sr = 5e6;
    double_t T = 1 / freq;
    float_t cycles = (measured_sec / T) / 2, data_size = sr * measured_sec;
    double_t x = 0, end_x = cycles * 2 * M_PI;

    double_t step_time = end_x / data_size, st = measured_sec / data_size;

    outputs->arr_size = (cycles * 2 * M_PI) / step_time;
    // std::cout << outputs->arr_size << std::endl;

    outputs->ys = (float_t *)malloc(outputs->arr_size * sizeof(float_t));
    outputs->ts = (float_t *)malloc(outputs->arr_size * sizeof(float_t));
    size_t i = 0;
    for (x; x <= end_x; x = x + step_time, i++)
    {
        *(outputs->ys + i) = f(x);
        outputs->average += *(outputs->ys + i);
        *(outputs->ts + i) = (st * i);
    }
    outputs->average /= outputs->arr_size;
    return outputs;
}
