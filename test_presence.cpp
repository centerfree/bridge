#include "uwb.h"
#include "stdio.h"

#define INPUT_FILE "/home/ashish/Desktop/uwb/x4dat/%d.bin"
#define OUTPUT_FILE "/home/ashish/Desktop/uwb/x4dat/%d.res"
#define MAX_DISTANCE 850
#define MIN_DISTANCE 200
#define INTERFERENCE_RANGE 167
#define BLOB_SIZE 40
#define SENSITIVITY 0.01

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
        fprintf(stdout, "%d, %d\n", height, width);

        float *array = (float *)malloc(sizeof(float) * width * height);
        read_float32_array(array, in_file);

        // Example of presence api call
        float result = presence(array, height, width, 100, MIN_DISTANCE, MAX_DISTANCE, INTERFERENCE_RANGE, BLOB_SIZE);

        printf("\nDetected presence with score %f", result);
    }
    return 0;
}
