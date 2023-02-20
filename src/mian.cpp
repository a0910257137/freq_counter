#include "../include/user_types.h"
#include "../include/utils.h"
#include "../include/zero_crossing.h"
#include "../include/interpolation.h"
#include "../include/fourier.h"
#include "../include/cmplx_funcs.h"
#include "matplotlibcpp.h"

using namespace std;
using std::vector;
namespace plt = matplotlibcpp;

int main(int argc, char *argv[])
{
    const unsigned int file_length = 250 * 1e6;
    const float Fs = 2.5 * 1e8;
    unsigned int glb_idx = 0, length = 96 * 96 * 96, rounds = 5;
    float *ys = (float *)malloc(length * sizeof(float)), *copy_ys = (float *)malloc(length * sizeof(float));
    // float f = 50 * 1e6;
    float f = 1000;
    FFT_data *fft_data = (FFT_data *)malloc(sizeof(FFT_data));
    zc_t *zc_data = (zc_t *)malloc(sizeof(zc_t));
    __init_fft_params(fft_data, length, Fs);
    __init_zc_params(zc_data, Fs);
    fft_data->Re = (float *)malloc(length * sizeof(float));
    while (rounds)
    {
        clock_t dft_start = clock();
        float val;
        for (int i = 0; i < length; i++)
        {
            val = sin(2 * PI * f * (i / (float)Fs));
            fft_data->Re[i] = val;
            copy_ys[i] = val;
        }
        improve_fft(fft_data->Re, fft_data->Im, fft_data->temp_Re, fft_data->temp_Im, fft_data->pre_build_table, fft_data->table_idx, length);
        cmplx_peak_mag(fft_data->Re, fft_data->Im, fft_data->magnitudes, fft_data->peak_magnitude, fft_data->peak,
                       length);
        fft_data->estimated_freq = fft_data->peak / fft_data->T;
        if (fft_data->estimated_freq < 1e5)
        {
            print(80);
            printf("Start applying zero-crossing method...\n");
            estimate(zc_data, copy_ys, length);
            if (zc_data->mean_freq == 0.)
                printf("The estimated frequency is  DC\n");
            else
                printf("The estimated frequency is  %5.5f kHz\n", zc_data->mean_freq / 1e3);
            __reset_zc_params(zc_data, Fs);
        }
        else
        {
            // apply parabolic or gaussian to fine-tune peak
            find_k_peaks(fft_data, 200);
            print(80);
            if ((fft_data->estimated_freq < 50 * 1e6))
                parabolic(fft_data, 9);

            else if (fft_data->estimated_freq > 50 * 1e6)
                gaussian(fft_data, 9);
            printf("The estimated frequency is %5.5f MHz\n", fft_data->estimated_freq / 1e6);
        }

        clock_t dft_end = clock();
        cout << "Elapsed time time: " << (dft_end - dft_start) * 1000000 / CLOCKS_PER_SEC << " us.\n";
        __reset_params(fft_data, length);
        rounds--;
    }
    __free_params(fft_data);
}