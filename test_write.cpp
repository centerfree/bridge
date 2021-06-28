#include "uwb.h"

#define OUTPUT_FILE "/home/ashish/Desktop/uwb/x4dat/test.bin"

int main(int argc, char const *argv[])
{
    float *array;
    int height, width;

    height = 3;
    width = 2;
    array = (float *)malloc(height * width * sizeof(float));

    for (size_t i = 0; i < height; i++)
    {
        for (size_t j = 0; j < width; j++)
        {
            array[i * width + j] = (float)(10 * i + j);
        }
    }

    write_float32_array(array, height, width, OUTPUT_FILE);

    return 0;
}
