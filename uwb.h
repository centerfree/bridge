#include <complex>

#define Complex             std::complex<float>


#define FILE_READ_ERROR     2
#define FILE_WRITE_ERROR    3
#define SUCCESS             1

#define DEFAULT_FC          (float)7.29e9
#define DEFAULT_FS          (float)23.328e9

// Array indexing: arr2d[row][col] = arr1d[row * n_col + col]

float max(float a, float b);

double max(double a, double b);

double mag(Complex a);

double mag(double real, double imag);

int lpf_zero_phase(float* signal_in, int size, float* signal_out, float ripple_db, float cutoff_freq_norm);

int lpf_zero_phase(double* signal_in, int size, double* signal_out, float ripple_db, float cutoff_freq_norm);

int lpf_zero_phase(Complex* signal_in, int size, Complex* signal_out, float ripple_db, float cutoff_freq_norm);

int read_float32_array(float *buffer, const char *filepath);

int read_float32_array_dims(int &time, int &distance, const char *filepath);

int write_float32_array(float *buffer, int time, int distance, const char *filepath);

int write_complex_array( Complex *buffer, int time, int distance, const char *filepath);

int demodulate(float *src, int time, int distance,  Complex *dst, float carrier_freq, float sampling_freq);

int demodulate(float *src, int time, int distance, float *dst, float carrier_freq, float sampling_freq);

int background_subtraction( Complex *src, int time, int distance,  Complex *dst);

int background_subtraction(float *src, int time, int distance, float *dst);

int decimate(float *big, int big_size, float *small, int factor);

int decimate( Complex *big, int big_size,  Complex *small, int factor);

float presence(float *src, int time, int distance, float fps,
               int min_distance_index, int max_distance_index, int interference_range_index, int blob_size);
