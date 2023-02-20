#include "../include/zero_crossing.h"
int sgn(float val)
{
    return (0 < val) - (val < 0);
}

int zero_crossing(float *y, int N)
{
    for (int i = 0, prev = sgn(y[0]); i <= N; i++)
    {
        int curr = sgn(y[i]);
        if (curr)
        {
            if (sgn(prev * curr) == -1)
            {
                return i;
            }
            prev = curr;
        }
    }
    return 0;
}

float refinement(int i, float *y, float time_step, unsigned int K, unsigned int inter_N)
{
    float tmp_ts[K], tmp_ys[K], inter_ts[inter_N], inter_ys[inter_N], output_time = 0;
    int j, m, rounds = 5, index = -1;
    for (j = 0; j < K; j++)
    {
        tmp_ts[j] = time_step * (i - K / 2 + j + 1);
        tmp_ys[j] = y[i - K / 2 + j + 1];
    }
    // interpolated points in sub time
    for (j = 0; j < rounds; j++)
    {
        float its = (tmp_ts[K - 1] - tmp_ts[0]) / inter_N, sub_time;
        for (m = 0; m <= inter_N; m++)
        {
            sub_time = tmp_ts[0] + (m * its);
            inter_ts[m] = sub_time;
            inter_ys[m] = Neville(tmp_ts, tmp_ys, sub_time, K);
            // print(80);
            // std::cout << inter_ts[m] << std::endl;
            // std::cout << inter_ys[m] << std::endl;
        }
        index = zero_crossing(inter_ys, inter_N);
        if (index == inter_N)
        {
            output_time = inter_ts[index];
            break;
        }

        else if (index == 0)
            break;

        else if ((index != 0) & (index != -1))
        {
            for (m = 0; m < K; m++)
            {
                tmp_ts[m] = inter_ts[(index - K / 2 + j + 1)];
                tmp_ys[m] = inter_ys[(index - K / 2 + j + 1)];
            }
        }
        else
        {
            output_time = inter_ts[index];
        }
    }
    // not found the best interpolated points
    if ((index == 0) & (output_time == 0.))
    {
        output_time = i * time_step;
    }

    return output_time;
    // // std::cout << i << std::endl;
    // std::cout << output_time << std::endl;
    // exit(1);
}

int get_avg_len(float *t, float *y, int N)
{
    float sampling_time = 4 * 1e-9;
    float *tmp1_ts = t, *tmp1_ys = y, prev_zc = 1, mark_start_time = 0, mark_end_time = 0;
    int zc_index = zero_crossing(tmp1_ys, N), avg_pts = 0; // estimate to get the first cycle

    float *tmp2_ts = tmp1_ts + zc_index, *tmp2_ys = tmp1_ys + zc_index;
    for (int i = 0; i < N; i++)
    {
        int zc_index = zero_crossing(tmp2_ys + i, N);

        if (i == 0)
            mark_start_time = *(tmp2_ts + i + zc_index);
        if (i != 0 & (zc_index / prev_zc) > 4)
        {
            mark_end_time = *(tmp2_ts + i + zc_index);
            float diffs = (mark_end_time - mark_start_time) / sampling_time;
            diffs < 0 ? -diffs : diffs;
            if (diffs > 300. & diffs < 350.)
                avg_pts = 35;
            else if (diffs > 250. & diffs < 300.)
                avg_pts = 30;
            else if (diffs > 200. & diffs < 250.)
                avg_pts = 25;
            else if (diffs > 150. & diffs < 200.)
                avg_pts = 20;
            else if (diffs > 100. & diffs < 150.)
                avg_pts = 15;

            else if (diffs > 50. & diffs < 100.)
                avg_pts = 10;

            else if (diffs > 10. & diffs < 50.)
                avg_pts = 5;
        }
        prev_zc = zc_index;
        if (avg_pts != 0)
            break;
    }

    return avg_pts;
}

float avg_val(float *val, int avg_len)
{
    float out_val = 0;
    for (int j = 0; j < avg_len; j++)
    {
        out_val += *(val + j);
    }

    return out_val / avg_len;
}
void estimate(zc_t *src_data, float *ys, const unsigned int &length)
{
    float curr_val, prev_val = *ys, prev_freq = 0, curr_freq = 0;
    int det_zc_times = 0, forward_step = 100;
    int prev = sgn(*ys), curr;
    for (int i = 0; i < length - 1;)
    {
        curr = sgn(ys[i]);
        if (curr)
        {
            curr_val = ys[i];
            if ((prev == +1) & (curr == -1) & sgn(prev * curr) == -1)
            {
                float interpolated_time = refinement(i, ys, src_data->time_step, src_data->K, src_data->inter_N);
                src_data->zc_ts[src_data->zc_counts] = interpolated_time;
                if (src_data->zc_counts != 0)
                {
                    det_zc_times++;
                    prev_freq = src_data->mean_freq;
                    curr_freq = 1 / ((src_data->zc_ts[src_data->zc_counts] - src_data->zc_ts[src_data->zc_counts - 1]));
                    src_data->freqs[det_zc_times - 1] = curr_freq;
                    src_data->mean_freq = prev_freq + (curr_freq - prev_freq) / det_zc_times;
                }
                if (src_data->zc_counts > buffer_size)
                    break;
                src_data->zc_counts++;
            }
            prev = curr;
            prev_val = curr_val;
        }

        i += forward_step;
    }
    // std::cout << src_data->mean_freq << std::endl;
    // exit(1);
}

void __init_zc_params(zc_t *src_data, const unsigned int &Fs)
{
    src_data->zc_counts = 0, src_data->mean_freq = 0, src_data->std_freq = 0;
    src_data->time_step = 1. / Fs;
    src_data->K = 5;
    src_data->inter_N = 20;
    src_data->zc_indices = (int *)malloc(buffer_size * sizeof(int));
    src_data->zc_ts = (float *)malloc(buffer_size * sizeof(float));
    src_data->freqs = (float *)malloc(buffer_size * sizeof(float));
}

void __reset_zc_params(zc_t *src_data, const unsigned int &Fs)
{
    src_data->zc_counts = 0, src_data->mean_freq = 0, src_data->std_freq = 0;
    src_data->time_step = 1. / Fs;
    src_data->K = 5;
    src_data->inter_N = 20;
    memset(src_data->zc_indices, 0., buffer_size * sizeof(float));
    memset(src_data->zc_ts, 0., buffer_size * sizeof(float));
    memset(src_data->freqs, 0., buffer_size * sizeof(float));
}