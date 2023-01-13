#include "../include/freq_counter.h"
namespace plt = matplotlibcpp;
template <typename T>
int sgn(T val)
{
    return (float(0) < val) - (val < float(0));
}

float Neville(float *ts, float *ys, float &t, size_t K)
{
    size_t i = 0;
    float val;

    float tmp_L[K] = {0}, tmp_R[K - 1] = {0};
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

float Linear(float &t0, float &t1, float &y0, float &y1, float &xt)
{

    return y0 + ((y1 - y0) / (t1 - t0)) * (xt - t0);
}
size_t zero_crossing(float_t *t, float_t *y, size_t &N)
{
    float_t curr_val, prev_val = *y;
    for (int i = 0, prev = sgn(*y); i < N - 1; i++)
    {
        if (int current = sgn(*(y + i)))
        {
            curr_val = *(y + i);
            if (sgn(prev * current) == -1)
            {
                return i;
            }
            prev = current;
        }
    }
    return 0;
}

void refinement(float_t *t, float_t *y, size_t i, size_t K, size_t inter_N)
{
    float_t time_range;
    float_t tmp_ts[K] = {0}, tmp_ys[K] = {0}, inter_ts[inter_N] = {0}, inter_ys[inter_N] = {0};

    for (int j = 0; j < K; j++)
    {
        tmp_ts[j] = *(t + (i - K / 2 + j + 1));
        tmp_ys[j] = *(y + (i - K / 2 + j + 1));
    }
    time_range = tmp_ts[K - 1] - tmp_ts[0];
    if (time_range != 0)
    {
        float its = time_range / inter_N;
        for (int l = 0; l <= inter_N; l++)
        {
            float_t t = tmp_ts[0] + l * its;
            inter_ts[l] = t;
            inter_ys[l] = Neville(tmp_ts, tmp_ys, t, K);
        }
        size_t index = zero_crossing(inter_ts, inter_ys, inter_N);
        //     float_t yp = Linear(t0, t1, y0, y1, xt);
        //     if (yp > 0)
        //     {
        //         xt = (yp - y0) * ((t1 - t0) / (y1 - y0)) + t0;
        //         y0 = yp;
        //         t0 = xt;
        //     }
        //     else
        //     {
        //         xt = (yp - y0) * ((t1 - t0) / (y1 - y0)) + t0;
        //         y1 = yp;
        //         t1 = xt;
        //     }
        //     xt = (t0 + t1) / 2;
        // }

        if (index != 0)
        {
            for (int j = 0; j < K; j++)
            {
                *(t + (i - K / 2 + j + 1)) = inter_ts[index - K / 2 + j + 1];
                *(y + (i - K / 2 + j + 1)) = inter_ys[index - K / 2 + j + 1];
            }
        }
    }
}

size_t random_number(int lower, int upper)
{
    return (rand() % (upper - lower + 1)) +
           lower;
}
size_t get_avg_len(float_t *t, float_t *y, size_t N)
{
    float_t sampling_time = 4 * 1e-9;
    float *tmp1_ts = t, *tmp1_ys = y, prev_zc = 1, mark_start_time = 0, mark_end_time = 0;
    size_t zc_index = zero_crossing(tmp1_ts, tmp1_ys, N), avg_pts; // estimate to get the first cycle
    float *tmp2_ts = tmp1_ts + zc_index, *tmp2_ys = tmp1_ys + zc_index;
    for (int i = 0; i < N; i++)
    {
        size_t zc_index = zero_crossing(tmp2_ts + i, tmp2_ys + i, N);
        if (i == 0)
            mark_start_time = *(tmp2_ts + i + zc_index);
        if (i != 0 & (zc_index / prev_zc) > 4)
        {
            mark_end_time = *(tmp2_ts + i + zc_index);
            float diffs = abs(mark_end_time - mark_start_time) / sampling_time;
            if (diffs > 300 & diffs < 350)
                avg_pts = 30;
            else if (diffs > 250 & diffs < 300)
                avg_pts = 25;
            else if (diffs > 200 & diffs < 250)
                avg_pts = 15;
            else if (diffs > 150 & diffs < 200)
                avg_pts = 20;
            else if (diffs > 100 & diffs < 150)
                avg_pts = 15;

            else if (diffs > 50 & diffs < 100)
                avg_pts = 10;

            else if (diffs > 10 & diffs < 50)
                avg_pts = 5;
        }
        prev_zc = zc_index;
        if (avg_pts != 0)
            break;
    }
    return avg_pts;
}

float_t avg_val(float_t *val, size_t &avg_len)
{
    float_t out_val = 0;
    for (int j = 0; j < avg_len; j++)
    {
        out_val += *(val + j);
    }

    return out_val / avg_len;
}
freq_t *estimate(point_t *data, size_t K, size_t inter_N)
{
    // select K around zero-crossing points and interpolated N points
    // sampling rate is 4ns for one point
    float_t curr_val, prev_val = *data->ys;
    freq_t *outputs = (freq_t *)malloc(1 * sizeof(freq_t));
    float_t prev_freq = 0, curr_freq = 0, mean_freq = 0, mean_T = 0, var_freq = 0, std_freq = 0, prec_in_t = 0;
    outputs->zc_counts = 0, outputs->mean_freq = 0, outputs->std_freq = 0, outputs->prec_in_t = 0;
    size_t freq_len = 0;
    size_t avg_pts = get_avg_len(data->ts, data->ys, data->length);
    // std::cout
    //     << avg_pts << std::endl;
    // exit(1);
    float_t *avg_ts = (float_t *)malloc(data->length * sizeof(float_t));
    float_t *avg_ys = (float_t *)malloc(data->length * sizeof(float_t));
    size_t avg_rounds = data->length;
    if (avg_pts != 0)
    {
        avg_rounds = data->length / avg_pts;

        for (int i = 0; i < avg_rounds; i++)
        {
            float_t avg_t = avg_val(data->ts + avg_pts * i, avg_pts);
            float_t avg_y = avg_val(data->ys + avg_pts * i, avg_pts);
            avg_ts[i] = avg_t, avg_ys[i] = avg_y;
        }
    }
    else
    {
        avg_ts = data->ts;
        avg_ys = data->ys;
    }
    for (int i = 0, prev = sgn(*avg_ys); i < avg_rounds - 1; i++)
    {
        if (int current = sgn(*(avg_ys + i)))
        {
            curr_val = *(avg_ys + i);
            if (prev == +1 & current == -1 & sgn(prev * current) == -1)
            {
                outputs->zc_indices[outputs->zc_counts] = i;

                // k times finetune
                for (int kk = 0; kk < 20; kk++)
                {
                    refinement(avg_ts, avg_ys, i, K, inter_N);
                }
                outputs->zc_ts[outputs->zc_counts] = (*(avg_ts + i - 1) + *(avg_ts + i)) / 2;
                outputs->zc_ys[outputs->zc_counts] = (*(avg_ys + i - 1) + *(avg_ys + i)) / 2;
                if (outputs->zc_counts != 0)
                {
                    freq_len++;

                    prev_freq = mean_freq;
                    curr_freq = 1 / ((outputs->zc_ts[outputs->zc_counts] - outputs->zc_ts[outputs->zc_counts - 1]));
                    outputs->freqs[freq_len - 1] = curr_freq;
                    mean_freq = prev_freq + (curr_freq - prev_freq) / freq_len;
                    mean_T = 1 / mean_freq;
                    prec_in_t = abs(mean_T * freq_len - outputs->zc_ts[outputs->zc_counts]);
                    if (freq_len != 1)
                        var_freq = (var_freq + abs(prev_freq - mean_freq) * abs(curr_freq - prev_freq));
                    // std_freq = sqrt(var_freq / freq_len);
                    // std::cout
                    //     << prec_in_t / freq_len << std::endl;
                }
                if (outputs->zc_counts > buffer_size)
                    break;
                outputs->zc_counts++;
            }
            prev = current;
            prev_val = curr_val;
        }
    }

    outputs->std_freq = sqrt(outputs->std_freq / freq_len);
    outputs->prec_in_t = prec_in_t / freq_len;
    outputs->mean_freq = mean_freq;
    // std::cout << mean_freq << std::endl;
    // std::cout << prec_in_t << std::endl;
    // exit(1);
    return outputs;
}
