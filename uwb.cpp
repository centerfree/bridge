#define _USE_MATH_DEFINES

#include "uwb.h"
#include "cmath"
#include "stdio.h"
#include "wavelib/header/wavelib.h"
#include "iir1/Iir.h"

float max(float a, float b)
{
    return a > b ? a : b;
}

double max(double a, double b)
{
    return a > b ? a : b;
}

double mag(Complex a)
{
    return sqrt(a.real() * a.real() + a.imag() * a.imag());
}

double mag(double real, double imag) { return sqrt(real * real + imag * imag); }

int read_float32_array(float *buffer, const char *filepath)
{
    FILE *array_file;
    int distance, time;

    if (!(array_file = fopen(filepath, "rt")))
    {
        printf("Error! opening file");
        return FILE_READ_ERROR;
    }

    float dims[2];
    fread(dims, sizeof(float), 2, array_file);
    time = (int)dims[0];
    distance = (int)dims[1];
    size_t count = fread(buffer, sizeof(float), distance * time, array_file);

    printf("%d read\n", (int)count);
    return SUCCESS;
}

int read_float32_array_dims(int &time, int &distance, const char *filepath)
{
    FILE *array_file;

    if (!(array_file = fopen(filepath, "rt")))
    {
        printf("Error! opening file");
        time = 0;
        distance = 0;
        return FILE_READ_ERROR;
    }

    float dims[2];
    fread(dims, sizeof(float), 2, array_file);
    time = (int)dims[0];
    distance = (int)dims[1];
    return SUCCESS;
}

int write_float32_array(float *buffer, int time, int distance,
                        const char *filepath)
{
    FILE *array_file;
    array_file = fopen(filepath, "wt");
    if (!array_file)
    {
        printf("Error! opening file.\n");
        return FILE_WRITE_ERROR;
    }

    float dims[2] = {(float)time, (float)distance};

    fwrite(dims, sizeof(float), 2, array_file);
    fwrite(buffer, sizeof(float), time * distance, array_file);

    // printf("Successfully written %d bytes.\n", 4 * time * distance);

    return SUCCESS;
}

int write_complex_array(Complex *buffer, int time, int distance,
                        const char *filepath)
{
    FILE *array_file;

    if (!(array_file = fopen(filepath, "wt")))
    {
        printf("Error! opening file.\n");
        return FILE_WRITE_ERROR;
    }

    float dims[2] = {(float)time, (float)distance};

    fwrite(dims, sizeof(float), 2, array_file);
    fwrite(buffer, sizeof(Complex), time * distance, array_file);

    return SUCCESS;
}

int demodulate(float *src, int time, int distance, Complex *dst,
               float carrier_freq, float sampling_freq)
{
    Complex *carrier_wave = (Complex *)malloc(distance * sizeof(Complex));
    Complex omega(0, -1);
    omega *= 2 * (float)M_PI * carrier_freq / sampling_freq;
    for (size_t i = 0; i < distance; i++)
    {
        // carrier = numpy.exp(-1j * 2 * numpy.pi * numpy.arange(arr.shape[1]) *
        // __FC / __FS)
        carrier_wave[i] = exp(omega * (Complex)(float)i);
    }

    for (size_t i = 0; i < time; i++)
    {
        for (size_t j = 0; j < distance; j++)
        {
            dst[i * distance + j] = (Complex)src[i * distance + j] * carrier_wave[j];
        }
    }
    free(carrier_wave);
    return SUCCESS;
}

int demodulate(float *src, int time, int distance, float *dst,
               float carrier_freq, float sampling_freq)
{
    Complex *carrier_wave = (Complex *)malloc(distance * sizeof(Complex));
    Complex omega(0, -1);
    omega *= 2 * (float)M_PI * carrier_freq / sampling_freq;
    for (size_t i = 0; i < distance; i++)
    {
        // carrier = numpy.exp(-1j * 2 * numpy.pi * numpy.arange(arr.shape[1]) *
        // __FC / __FS)
        carrier_wave[i] = exp(omega * (Complex)(float)i);
    }

    for (size_t i = 0; i < time; i++)
    {
        for (size_t j = 0; j < distance; j++)
        {
            dst[i * distance + j] =
                (float)mag(((Complex)src[i * distance + j] * carrier_wave[j]));
        }
    }
    free(carrier_wave);
    return SUCCESS;
}

int background_subtraction(Complex *src, int time, int distance,
                           Complex *dst)
{
    Complex background;

    for (size_t d = 0; d < distance; d++)
    {
        background = 0;
        for (size_t t = 0; t < time; t++)
        {
            background += src[t * distance + d];
        }
        background /= (float)time;

        for (size_t t = 0; t < time; t++)
        {
            dst[t * distance + d] = src[t * distance + d] - background;
        }
    }
    return SUCCESS;
}

int background_subtraction(float *src, int time, int distance, float *dst)
{
    double background;

    for (size_t d = 0; d < distance; d++)
    {
        background = 0;
        for (size_t t = 0; t < time; t++)
        {
            background += (double)src[t * distance + d];
        }
        background /= (double)time;

        for (size_t t = 0; t < time; t++)
        {
            dst[t * distance + d] = src[t * distance + d] - (float)background;
        }
    }
    return SUCCESS;
}


int lpf_zero_phase(float *signal_in, int size, float *signal_out, float ripple_db, float cutoff_freq_norm)
{
    const int order = 8;
    Iir::ChebyshevI::LowPass<order> lpf;
    lpf.setupN(0.5 * cutoff_freq_norm, ripple_db);
    Iir::ChebyshevI::LowPass<order> lpb;
    lpb.setupN(0.5 * cutoff_freq_norm, ripple_db);

    int pad = (int)(size * cutoff_freq_norm);

    float *tmp = (float *)malloc(sizeof(float) * (size + 2 * pad));

    for (size_t i = 0; i < pad; i++)
        tmp[i] = lpf.filter(signal_in[pad - i]);

    for (size_t i = 0; i < size; i++)
        tmp[i + pad] = lpf.filter(signal_in[i]);

    for (size_t i = 0; i < pad; i++)
        tmp[i + pad + size] = lpf.filter(signal_in[size - i - 2]);

    for (size_t i = 0; i < pad; i++)
        signal_out[0] = lpb.filter(tmp[size + 2 * pad - i - 1]);

    for (size_t i = 0; i < size; i++)
        signal_out[size - i - 1] = lpb.filter(tmp[size + pad - 1 - i]);

    free(tmp);
    return SUCCESS;
}

int lpf_zero_phase(double *signal_in, int size, double *signal_out, float ripple_db, float cutoff_freq_norm)
{
    const int order = 8;
    Iir::ChebyshevI::LowPass<order> lpf;
    lpf.setupN(0.5 * cutoff_freq_norm, ripple_db);
    Iir::ChebyshevI::LowPass<order> lpb;
    lpb.setupN(0.5 * cutoff_freq_norm, ripple_db);

    int pad = (int)(size * cutoff_freq_norm);

    float *tmp = (float *)malloc(sizeof(float) * (size + 2 * pad));

    for (size_t i = 0; i < pad; i++)
        tmp[i] = lpf.filter((float)signal_in[pad - i]);

    for (size_t i = 0; i < size; i++)
        tmp[i + pad] = lpf.filter((float)signal_in[i]);

    for (size_t i = 0; i < pad; i++)
        tmp[i + pad + size] = lpf.filter((float)signal_in[size - i - 2]);

    for (size_t i = 0; i < pad; i++)
        signal_out[0] = (double)lpb.filter(tmp[size + 2 * pad - i - 1]);

    for (size_t i = 0; i < size; i++)
        signal_out[size - i - 1] = (double)lpb.filter(tmp[size + pad - 1 - i]);

    free(tmp);
    return SUCCESS;
}

int lpf_zero_phase(Complex *signal_in, int size, Complex *signal_out, float ripple_db, float cutoff_freq_norm)
{
    const int order = 8; // 4th order (=2 biquads)
    Iir::ChebyshevI::LowPass<order> realf;
    realf.setupN(0.5 * cutoff_freq_norm, ripple_db);
    Iir::ChebyshevI::LowPass<order> realb;
    realb.setupN(0.5 * cutoff_freq_norm, ripple_db);
    Iir::ChebyshevI::LowPass<order> imagf;
    imagf.setupN(0.5 * cutoff_freq_norm, ripple_db);
    Iir::ChebyshevI::LowPass<order> imagb;
    imagb.setupN(0.5 * cutoff_freq_norm, ripple_db);

    int pad = (int)(size * cutoff_freq_norm);

    Complex *tmp = (Complex *)malloc(sizeof(Complex) * (size + 2 * pad));

    for (size_t i = 0; i < pad; i++)
        tmp[i] = Complex(realf.filter(signal_in[pad - i].real()), imagf.filter(signal_in[pad - i].imag()));

    for (size_t i = 0; i < size; i++)
        tmp[i + pad] = Complex(realf.filter(signal_in[i].real()), imagf.filter(signal_in[i].imag()));

    for (size_t i = 0; i < pad; i++)
        tmp[i + pad + size] = Complex(realf.filter(signal_in[size - i - 2].real()), imagf.filter(signal_in[size - i - 2].imag()));

    for (size_t i = 0; i < pad; i++)
        signal_out[0] = Complex(realb.filter(tmp[size + 2 * pad - i - 1].real()), imagb.filter(tmp[size + 2 * pad - i - 1].imag()));

    for (size_t i = 0; i < size; i++)
        signal_out[size - i - 1] = Complex(realb.filter(tmp[size + pad - 1 - i].real()), imagb.filter(tmp[size + pad - 1 - i].imag()));

    free(tmp);
    return SUCCESS;
}

int decimate(float *big, int big_size, float *small, int factor)
{
    float *trim = (float *)malloc(sizeof(float) * (big_size));
    lpf_zero_phase(big, big_size, trim, 0.05, 0.8 / (float)factor);

    int j = 0;
    while (j * factor < big_size)
    {
        small[j] = trim[j * factor];
        j += 1;
    }

    free(trim);
    return SUCCESS;
}

int decimate(Complex *big, int big_size, Complex *small, int factor)
{
    Complex *trim = (Complex *)malloc(sizeof(Complex) * (big_size));
    lpf_zero_phase(big, big_size, trim, 0.05, 0.8 / (float)factor);

    int j = 0;
    while (j * factor < big_size)
    {
        small[j] = trim[j * factor];
        j += 1;
    }

    free(trim);
    return SUCCESS;
}

float presence(float *src, int time, int distance, float fps,
               int min_distance_index, int max_distance_index,
               int interference_range_index, int blob_size)
{

    Complex *processed =
        (Complex *)malloc(time * distance * sizeof(Complex));
    Complex *clean = (Complex *)malloc(time * distance * sizeof(Complex));

    demodulate(src, time, distance, processed, DEFAULT_FC, DEFAULT_FS);
    background_subtraction(processed, time, distance, clean);

    char wavelet_name[16] = "db2";
    wave_object wavelet = wave_init(wavelet_name);

    double *signal_real = (double *)malloc(sizeof(double) * time);
    double *signal_imag = (double *)malloc(sizeof(double) * time);

    float *raw_energy = (float *)malloc(
        sizeof(float) * (max_distance_index - min_distance_index + 1));

    int max_levels = (int)log2((double)time / (double)wavelet->filtlength);

    wt_object wt_real = wt_init(wavelet, "dwt", time, max_levels);
    wt_object wt_imag = wt_init(wavelet, "dwt", time, max_levels);

    // printf("%d filter coeffs, %d levels for array sized %d\n",
    // wavelet->filtlength, max_levels, time);

    for (int d = min_distance_index; d <= max_distance_index; d++)
    {
        for (int t = 0; t < time; t++)
        {
            signal_real[t] = clean[t * distance + d].real();
            signal_imag[t] = clean[t * distance + d].imag();
        }

        dwt(wt_real, signal_real);
        dwt(wt_imag, signal_imag);

        raw_energy[d - min_distance_index] = 0;
        for (size_t c = 0; c < wt_real->length[0]; c++)
        {
            raw_energy[d - min_distance_index] +=
                (float)mag(wt_real->output[c], wt_imag->output[c]);
        }
        raw_energy[d - min_distance_index] /= wt_real->length[0];
    }

    float e, max_e;
    int left, right;
    max_e = 0;
    for (int d = min_distance_index; d <= max_distance_index; d++)
    {
        e = raw_energy[d - min_distance_index];
        for (int i = 1; i < int(blob_size / 2); i++)
        {
            left = d - i;
            right = d + i;
            left = left >= min_distance_index ? left : 2 * min_distance_index - left;
            right =
                right <= max_distance_index ? right : 2 * max_distance_index - right;
            e += raw_energy[left - min_distance_index] +
                 raw_energy[right - min_distance_index];
        }
        e /= (float)(blob_size + 1);
        max_e = max(max_e, e);
    }

    free(processed);
    free(clean);
    free(signal_imag);
    free(signal_real);
    free(raw_energy);

    return (float)max_e;
}

