#include "iir1/Iir.h"
#include "uwb.h"
#include "math.h"
#include "malloc.h"

int main(int argc, char const *argv[])
{
    float start = 0;
    float rate = 0.5;
    int len = 1000;
    int factor = 4;

    float *signal = (float *)malloc(sizeof(float) * len);
    float *decimated = (float *)malloc(sizeof(float) * int(len / factor));

    for (size_t i = 0; i < len; i++)
    {
        signal[i] = cos((float)(i * (start + i * rate)) / (float)len);
    }

    decimate(signal, len, decimated, factor);

    write_float32_array(signal, 1, len, "../x4dat/rc1.txt");
    write_float32_array(decimated, 1, int(len / factor), "../x4dat/rcd1.txt");

    Complex *signal_complex = (Complex *)malloc(sizeof(Complex) * len);
    Complex *decimated_complex = (Complex *)malloc(sizeof(Complex) * int(len / factor));

    for (size_t i = 0; i < len; i++)
    {
        signal_complex[i] = exp( Complex(0, 1) * (float)(i * (start + i * rate) / len));
    }

    decimate(signal_complex, len, decimated_complex, factor);

    write_complex_array(signal_complex, 1, len, "../x4dat/cc1.txt");
    write_complex_array(decimated_complex, 1, int(len / factor), "../x4dat/ccd1.txt");

    return 0;
}