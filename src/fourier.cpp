#include "../include/fourier.h"

void seperate(float *Re, float *Im, float *ve_Re, float *ve_Im, float *vo_Re, float *vo_Im, unsigned int n)
{
    int blkCnt = n >> 2; // 1/2 for even and odd; the other is wating calculation
    for (int k = 0; k < blkCnt; k++)
    {
        *ve_Re++ = *Re;
        *ve_Im++ = *Im;
        Re++;
        Im++;
        *vo_Re++ = *Re;
        *vo_Im++ = *Im;
        Re++;
        Im++;

        *ve_Re++ = *Re;
        *ve_Im++ = *Im;

        Re++;
        Im++;

        *vo_Re++ = *Re;
        *vo_Im++ = *Im;
        Re++;
        Im++;
    }
    // *ve_Re++ = *Re;
    // *ve_Im++ = *Im;
    // Re++;
    // Im++;
    // *vo_Re++ = *Re;
    // *vo_Im++ = *Im;
    // Re++;
    // Im++;
}

void _get_common_table(float *table, unsigned int &glb_idx, unsigned int n)
{
    if (n > 1)
    {
        int m = 0;
        int l = n >> 1;
        _get_common_table(table, glb_idx, l);
        _get_common_table(table, glb_idx, l);
        float val_cos, val_sin;
        for (m = 0; m < l; m++)
        {
            *(table + glb_idx) = cos(2 * PI * m / n);
            glb_idx += 1;
            *(table + glb_idx) = -sin(2 * PI * m / n);
            glb_idx += 1;
        }
    }
    return;
}
void _calculate(float *src_Re, float *src_Im, float *vo_Re, float *vo_Im, float *ve_Re, float *ve_Im, float *table, unsigned int &glb_idx, unsigned int &l)
{
    float z_Re, z_Im, w_Re, w_Im;
    float vo1 = 0, vo2 = 0, ve1 = 0, ve2 = 0;
    for (int m = 0; m < l; m++)
    {
        w_Re = table[glb_idx];
        glb_idx += 1;
        w_Im = table[glb_idx];
        glb_idx += 1;
        vo1 = *vo_Re++;
        vo2 = *vo_Im++;
        ve1 = *ve_Re++;
        ve2 = *ve_Im++;
        z_Re = w_Re * vo1 - w_Im * vo2; /* Re(w*vo[m]) */
        z_Im = w_Re * vo2 + w_Im * vo1; /* Im(w*vo[m]) */
        *src_Re++ = ve1 + z_Re;
        *src_Im++ = ve2 + z_Im;
        src_Re[l - 1] = ve1 - z_Re;
        src_Im[l - 1] = ve2 - z_Im;
    }
}
void improve_fft(float *src_Re, float *src_Im, float *temp_Re, float *temp_Im, float *table, unsigned int &glb_idx, unsigned int n)
{
    if (n > 1)
    {
        unsigned int l, m;
        l = n >> 1;
        float *vo_Re, *vo_Im, *ve_Re, *ve_Im;
        ve_Re = temp_Re;
        ve_Im = temp_Im;
        vo_Re = temp_Re + l;
        vo_Im = temp_Im + l;
        seperate(src_Re, src_Im, ve_Re, ve_Im, vo_Re, vo_Im, n);
        improve_fft(ve_Re, ve_Im, src_Re, src_Im, table, glb_idx, l); /* FFT on even-indexed elements of v[] */
        improve_fft(vo_Re, vo_Im, src_Re, src_Im, table, glb_idx, l); /* FFT on odd-indexed elements of v[] */
        _calculate(src_Re, src_Im, vo_Re, vo_Im, ve_Re, ve_Im, table, glb_idx, l);
    }
    return;
}

void fft(complex *v, int n, complex *tmp)
{
    if (n > 1)
    {
        int k, m;
        complex z, w, *vo, *ve;
        ve = tmp;
        vo = tmp + n / 2;
        for (k = 0; k < n / 2; k++)
        {
            *(ve + k) = *(v + 2 * k);
            *(vo + k) = *(v + 2 * k + 1);
        }
        fft(ve, n / 2, v); /* FFT on even-indexed elements of v[] */
        fft(vo, n / 2, v); /* FFT on odd-indexed elements of v[] */
        for (m = 0; m < n / 2; m++)
        {
            w.Re = cos(2 * PI * m / n);
            w.Im = -sin(2 * PI * m / n);
            z.Re = w.Re * vo[m].Re - w.Im * vo[m].Im; /* Re(w*vo[m]) */
            z.Im = w.Re * vo[m].Im + w.Im * vo[m].Re; /* Im(w*vo[m]) */
            v[m].Re = ve[m].Re + z.Re;
            v[m].Im = ve[m].Im + z.Im;
            v[m + n / 2].Re = ve[m].Re - z.Re;
            v[m + n / 2].Im = ve[m].Im - z.Im;
        }
    }
    return;
}

void __init_fft_params(FFT_data *src_data, const unsigned int &length, const unsigned int &Fs)
{
    src_data->top_k = 9;
    src_data->pre_build_table = (float *)malloc(length * 21 * sizeof(float));
    src_data->Re = (float *)malloc(length * sizeof(float));
    src_data->Im = (float *)calloc(length, sizeof(float));
    src_data->temp_Re = (float *)malloc(length * sizeof(float));
    src_data->temp_Im = (float *)malloc(length * sizeof(float));
    src_data->magnitudes = (float *)malloc(length * sizeof(float));
    src_data->T = (float)length / (float)Fs;
    src_data->k_frequencies = (float *)malloc(src_data->top_k * sizeof(float));
    src_data->k_maginitudes = (float *)malloc(src_data->top_k * sizeof(float));
    src_data->k_peaks = (float *)malloc(src_data->top_k * sizeof(float));
    src_data->peak = 0;
    src_data->table_idx = 0;
    _get_common_table(src_data->pre_build_table, src_data->table_idx, length);
    src_data->table_idx = 0;
}

void __reset_params(FFT_data *src_data, const unsigned int &length)
{
    memset(src_data->Re, 0., length * sizeof(float));
    memset(src_data->Im, 0., length * sizeof(float));
    memset(src_data->magnitudes, 0., length * sizeof(float));
    memset(src_data->temp_Re, 0., length * sizeof(float));
    memset(src_data->temp_Im, 0., length * sizeof(float));
    memset(src_data->k_frequencies, 0., src_data->top_k * sizeof(float));
    memset(src_data->k_maginitudes, 0., src_data->top_k * sizeof(float));
    memset(src_data->k_peaks, 0., src_data->top_k * sizeof(float));
    src_data->peak = 0;
    src_data->table_idx = 0;
}

void __free_params(FFT_data *src_data)
{
    free(src_data->Re);
    free(src_data->Im);
    free(src_data->magnitudes);
    free(src_data->temp_Re);
    free(src_data->temp_Im);
    free(src_data->k_frequencies);
    free(src_data->k_maginitudes);
    free(src_data->k_peaks);
}