#include "../include/cmplx_funcs.h"

void cmplx_peak_mag(float *Re,
                    float *Im,
                    float *pDst,
                    float &magnitude,
                    unsigned int &peak,
                    const unsigned int numSamples)
{

    unsigned int blkCnt, final_peak, i = 0, tol = numSamples >> 1;
    float realIn, imagIn, max_val, final_max_val; /* Temporary variables to hold input values */
    blkCnt = numSamples >> 2u;
    while (blkCnt > 0u)
    {
        /* C[0] = sqrt(A[0] * A[0] + A[1] * A[1]) */
        realIn = *Re++;
        imagIn = *Im++;
        pDst++;
        /* store the result in the destination buffer. */
        *pDst = sqrt((realIn * realIn) + (imagIn * imagIn));
        if ((*pDst > max_val) & (i != 0) & (i != numSamples - 1))
        {
            max_val = *pDst;
            peak = i;
        }
        i++;

        realIn = *Re++;
        imagIn = *Im++;
        pDst++;
        *pDst = sqrt((realIn * realIn) + (imagIn * imagIn));
        if ((*pDst > max_val) & (i != 0) & (i != numSamples - 1))
        {
            max_val = *pDst;
            peak = i;
        }
        i++;
        realIn = *Re++;
        imagIn = *Im++;
        pDst++;
        *pDst = sqrt((realIn * realIn) + (imagIn * imagIn));
        if ((*pDst > max_val) & (i != 0) & (i != numSamples - 1))
        {
            max_val = *pDst;
            peak = i;
        }
        i++;
        realIn = *Re++;
        imagIn = *Im++;
        pDst++;
        *pDst = sqrt((realIn * realIn) + (imagIn * imagIn));
        if ((*pDst > max_val) & (i != 0) & (i != numSamples - 1))
        {
            max_val = *pDst;
            peak = i;
        }
        i++;
        /* Decrement the loop counter */
        blkCnt--;
        if ((i <= numSamples >> 1 + tol)) // to know where is maximum
        {
            final_max_val = max_val;
            final_peak = peak;
            // break;
        }
    }
    magnitude = final_max_val;
    peak = final_peak;
}

void find_k_peaks(FFT_data *src,
                  const unsigned int searching_range)
{
    unsigned int i, j, top_k = 9, L = top_k / 2, l = 0;
    float *buffer = (float *)malloc(top_k * sizeof(float)), value;
    src->k_peaks[L] = src->peak;
    src->k_frequencies[L] = (float)src->peak / src->T;
    src->k_maginitudes[L] = src->magnitudes[src->peak];
    // LHS
    for (i = 1; i < searching_range >> 1; i++)
    {
        value = src->magnitudes[src->peak - i];
        for (j = L; j != 0; j--)
        {
            if (src->k_maginitudes[j - 1] < value)
            {
                src->k_peaks[j - 1] = src->peak - i;
                src->k_frequencies[j - 1] = (float)(src->peak - i) / src->T;
                src->k_maginitudes[j - 1] = value;
                break;
            }
        }
    }
    // RHS
    for (i = 1; i < searching_range / 2; i++)
    {
        value = src->magnitudes[src->peak + i];
        for (j = L; j < top_k; j++)
        {
            if (src->k_maginitudes[j + 1] < value)
            {
                src->k_peaks[j + 1] = src->peak + i;
                src->k_frequencies[j + 1] = (float)(src->peak + i) / src->T;
                src->k_maginitudes[j + 1] = value;
                break;
            }
        }
    }
}
