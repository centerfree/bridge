#include "uwb.h"
#include "stdio.h"

#define INPUT_FILE "/home/ashish/Desktop/uwb/x4dat/test_np.bin"

int main(int argc, char const *argv[])
{
    int height, width;

    read_float32_array_dims(height, width, INPUT_FILE);
    fprintf(stdout, "%d, %d\n\n", height, width);

    float *array = (float *)malloc(sizeof(float) * width * height);
    read_float32_array(array, INPUT_FILE);

    for (size_t i = 0; i < height; i++)
    {
        for (size_t j = 0; j < width; j++)
        {
            fprintf(stdout, "%d: %f\n", (int)(i * width + j), array[i * width + j]);
        }
    }

    return 0;
}
