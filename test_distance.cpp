#include "uwb.h"
#include "distance.h"
#include "stdio.h"
#include "malloc.h"

#define INPUT_FILE "../x4dat/%d.bin"
#define OUTPUT_FILE "../x4dat/%d.dis"
#define MAX_DISTANCE 850
#define MIN_DISTANCE 25
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
        // float *energy = (float *)malloc(sizeof(float) * (MAX_DISTANCE - MIN_DISTANCE + 1));
        float *energy = (float *)malloc(sizeof(float) * width);
        
        read_float32_array(array, in_file);

        distance(array, height, width, 100.0f, energy, MIN_DISTANCE, MAX_DISTANCE, INTERFERENCE_RANGE,
                 SENSITIVITY, BLOB_SIZE);

        write_float32_array(energy, 1, width, out_file);
    }
    return 0;
}
