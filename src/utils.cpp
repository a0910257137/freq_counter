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
std::vector<T> to_vector(size_t n, T(*x))
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

float_t *_read(std::string &path, size_t length, FILE *ptr)
{

    float_t *buffer = (float_t *)malloc(length * sizeof(float_t));
    char char_path[path.length() + 1];
    strcpy(char_path, path.c_str());
    ptr = fopen(char_path, "rb");
    fread(buffer, length * sizeof(float_t), 1, ptr);
    return buffer;
}

void do_fetch_data(raw_data_t *raw_data, point_t *points, const size_t &size)
{
    // unsigned short int j = abs(rand()) / 1e7;
    unsigned short int j = 0;

    std::memcpy(points->ts, (raw_data->ts + j), size * sizeof(float_t));
    std::memcpy(points->ys, (raw_data->ys + j), size * sizeof(float_t));

    // for (int i = 0; i < points->length; i++)
    // {
    //     *(points->ts + i) = *(raw_data->ts + i + j);
    //     *(points->ys + i) = *(raw_data->ys + i + j);
    // }

    // return points;
}