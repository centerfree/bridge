#include "uwb.h"
#include "stdio.h"

#define INPUT_FILE "/home/ashish/Desktop/uwb/x4dat/%d.bin"
#define OUTPUT_FILE "/home/ashish/Desktop/uwb/x4dat/%d.clean"

int main(int argc, char const *argv[])
{
    int height, width;

    for (size_t i = 1; i < 4; i++)
    {
        char in_file[256];
        char out_file[256];
        sprintf(in_file, INPUT_FILE, (int)i);
        sprintf(out_file, OUTPUT_FILE, (int)i);

        read_float32_array_dims(height, width, in_file);
        fprintf(stdout, "%d, %d\n\n", height, width);

        float *array = (float *)malloc(sizeof(float) * width * height);
        read_float32_array(array, in_file);

        float *processed = (float *)malloc(height * width * sizeof(float));

        demodulate(array, height, width, processed, DEFAULT_FC, DEFAULT_FS);
        background_subtraction(processed, height, width, array);

        write_float32_array(array, height, width, out_file);
    }

    return 0;
}
