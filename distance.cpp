#include "uwb.h"
#include "distance.h"
#include "wavelib/header/wavelib.h"

float distance(float *src, int time, int distance, float fps, float *energy,
               int min_distance_index, int max_distance_index,
               int interference_range_index, float sensitivity, int blob_size)
{
    Complex *processed = (Complex *)malloc(time * distance * sizeof(Complex));
    Complex *clean = (Complex *)malloc(time * distance * sizeof(Complex));

    demodulate(src, time, distance, processed, DEFAULT_FC, (float)DEFAULT_FS);
    background_subtraction(processed, time, distance, clean);

    char wavelet_name[16] = "db2";
    wave_object wavelet = wave_init(wavelet_name);

    double *signal_real = (double *)malloc(sizeof(double) * time);
    double *signal_imag = (double *)malloc(sizeof(double) * time);

    double *raw_energy = (double *)malloc(
        sizeof(double) * (max_distance_index - min_distance_index + 1));
    double *filtered_energy = (double *)malloc(
        sizeof(double) * (max_distance_index - min_distance_index + 1));

    int max_levels = (int)log2((double)time / (double)wavelet->filtlength);

    wt_object wt_real = wt_init(wavelet, "dwt", time, max_levels);
    wt_object wt_imag = wt_init(wavelet, "dwt", time, max_levels);

    // printf("%d filter coeffs, %d levels for array sized %d\n",
    // wavelet->filtlength, max_levels, time);

    for (size_t d = min_distance_index; d <= max_distance_index; d++)
    {
        for (size_t t = 0; t < time; t++)
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
                mag(wt_real->output[c], wt_imag->output[c]);
        }
        raw_energy[d - min_distance_index] /= wt_real->length[0];
    }

    lpf_zero_phase(raw_energy, max_distance_index - min_distance_index + 1, filtered_energy, 0.05, 0.2);

    // scale for adjusting rx-tx interference
    if (interference_range_index > min_distance_index)
    {
        for (size_t d = min_distance_index; d < interference_range_index; d++)
        {
            filtered_energy[d - min_distance_index] *= (double)(d * d) / (double)(interference_range_index * interference_range_index);
        }
    }

    for (size_t d = 0; d < min_distance_index; d++)
        energy[d] = 0;

    for (size_t d = max_distance_index + 1; d < distance; d++)
        energy[d] = 0;
    
    for (size_t d = min_distance_index; d <= max_distance_index; d++){
        energy[d] = (float)filtered_energy[d - min_distance_index];
    }

    double e;
    int index_e = -1;
    int left, right;
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
        e /= (double)(blob_size + 1);
        index_e = e > sensitivity ? d : index_e;
    }

    free(processed);
    free(clean);
    free(signal_imag);
    free(signal_real);
    free(raw_energy);

    return (float)index_e;
}
