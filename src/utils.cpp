#include "../include/utils.h"
void print(int len)
{
    for (int i = 0; i < len; i++)
    {
        printf("*");
    }
    printf("\n");
}
float diff_ms(struct timespec *t1, struct timespec *t0)
{
    struct timespec d = {.tv_sec = t1->tv_sec - t0->tv_sec, .tv_nsec = t1->tv_nsec - t0->tv_nsec};
    if (d.tv_nsec < 0)
    {
        d.tv_nsec += 1000000000;
        d.tv_sec--;
    }
    return (double)(d.tv_nsec / 1000000);
}
template <typename T>
std::vector<T> to_vector(int n, T(*x))
{
    std::vector<T> output_vec(n, 0);
    for (int i = 0; i < n; i++)
    {
        output_vec[i] = *(x + i);
    }
    return output_vec;
}
std::vector<double> to_vectorf64(int n, double(*x))
{
    return to_vector(n, x);
}
std::vector<float> to_vectorf32(int n, float(*x))
{
    return to_vector(n, x);
}

float *_read(const char *path, int length, FILE *ptr)
{

    float *buffer = (float *)malloc(length * sizeof(float));
    ptr = fopen(path, "rb");
    fread(buffer, length * sizeof(float), 1, ptr);
    return buffer;
}

void do_fetch_data(float *raw_data, float *ys, unsigned int size)
{
    memcpy(ys, raw_data, size * sizeof(float));
}

void min_max_norm(float *ys, int N)
{
    float min_val = 0, max_val = 0, value;
    unsigned int i;
    for (i = 0; i < N; i++)
    {
        value = *ys;
        if (*ys > max_val)
            max_val = value;
        if (*ys < min_val)
            min_val = value;
        ys++;
    }
    for (i = 0; i < N; i++)
    {
        *ys = (*ys - min_val) / (max_val - min_val);
        ys--;
    }
}