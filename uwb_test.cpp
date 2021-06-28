#include "uwb.h"

#define INPUT_FILE  "/home/ashish/Desktop/1.bin"
#define OUTPUT_FILE "/home/ashish/Desktop/2.bin"

int main(int argc, char const *argv[])
{
    float* array;
    int height, width;
    read_float32_array(array, height, width, INPUT_FILE);

    float* processed = (float*) malloc(height * width * sizeof(float));

    demodulate(array, height, width, processed, DEFAULT_FC, DEFAULT_FS);

    write_float32_array(processed, height, width, OUTPUT_FILE);

    return 0;
}
