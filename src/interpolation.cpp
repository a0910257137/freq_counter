#include "../include/interpolation.h"

void parabolic(FFT_data *src, int top_k)
{
    int rounds = 1;
    int idx = top_k / 2;
    float ym1 = src->k_maginitudes[idx - 1];
    float y0 = src->k_maginitudes[idx];
    float ypl = src->k_maginitudes[idx + 1];
    float p = 0.5 * (ypl - ym1) / ((2 * y0 - (ypl + ym1)));
    float y = y0 - 0.25 * (ym1 - ypl) * p;
    src->estimated_freq = src->k_frequencies[4] + (p / src->T);
}

void gaussian(FFT_data *src, int top_k)
{
    int idx = top_k / 2;
    float ym1, y0, ypl;
    ym1 = src->k_maginitudes[idx - 1];
    y0 = src->k_maginitudes[idx];
    ypl = src->k_maginitudes[idx + 1];
    float p = 0.5 * (log(ypl) - log(ym1)) / (0.5 * (2 * log(y0) - log(ypl) - log(ym1)));
    // std::cout << (p / src->T) << std::endl;
    // exit(1);
    src->estimated_freq = src->k_frequencies[4] + (p / src->T);
}
float Neville(float *ts, float *ys, float t, int K)
{
    int i = 0;
    float val;
    float tmp_L[K], tmp_R[K - 1];
    for (int k = 0; k < K; k++)
    {
        tmp_L[k] = *(ys + k);
    }

    while (K > 1)
    {
        for (int k = 0; k < K - 1; k++)
        {
            val = ((t - ts[k]) * tmp_L[k + 1] - (t - ts[k + 1 + i]) * tmp_L[k]) / (ts[k + 1 + i] - ts[k]);
            tmp_R[k] = val;
        }
        for (int k = 0; k < K - 1; k++)
            tmp_L[k] = tmp_R[k];
        i++;
        K--;
    }
    return tmp_L[0];
}