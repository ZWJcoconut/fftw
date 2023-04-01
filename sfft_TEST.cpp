#include <vector>
#include <complex>
#include <cmath>
#include <fftw3.h>
using namespace std;



std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& input) {
    int N = input.size();
    fftw_complex *in = reinterpret_cast<fftw_complex*>(fftw_alloc_complex(N));
    fftw_complex *out = reinterpret_cast<fftw_complex*>(fftw_alloc_complex(N));

    for (int i = 0; i < N; i++) {
        in[i][0] = input[i].real();
        in[i][1] = input[i].imag();
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    std::vector<std::complex<double>> output(N);
    for (int i = 0; i < N; i++) {
        output[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return output;
}

std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& input) {
    int N = input.size();
    fftw_complex *in = reinterpret_cast<fftw_complex*>(fftw_alloc_complex(N));
    fftw_complex *out = reinterpret_cast<fftw_complex*>(fftw_alloc_complex(N));

    for (int i = 0; i < N; i++) {
        in[i][0] = input[i].real();
        in[i][1] = input[i].imag();
    }

    fftw_plan plan = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    std::vector<std::complex<double>> output(N);
    for (int i = 0; i < N; i++) {
        output[i] = std::complex<double>(out[i][0], out[i][1]) / static_cast<double>(N);
    }

    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);

    return output;
}

std::vector<std::vector<std::complex<double>>> sfft(const std::vector<std::vector<std::complex<double>>>& X_nm) {
    int N = X_nm.size();
    int M = X_nm[0].size();

    std::vector<std::vector<std::complex<double>>> x_kl(N, std::vector<std::complex<double>>(M));
    std::vector<std::complex<double>> temp(M);

    double scale = sqrt(1.0 / N / M) * N;

    for (int k = 0; k < N; k++) {
        for (int l = 0; l < M; l++) {
            temp[l] = X_nm[k][l];
        }
        std::vector<std::complex<double>> fft_temp = fft(temp);
        std::vector<std::complex<double>> ifft_temp = ifft(fft_temp);
        for (int l = 0; l < M; l++) {
            x_kl[k][l] = scale * ifft_temp[l];
        }
    }

    return x_kl;
}
