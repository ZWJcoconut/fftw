#include <tuple>
#include <random>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <iostream>
#include <fftw3.h>
#include <cassert>


std::vector<std::complex<double>> fft(const std::vector<std::complex<double>>& vec) {
    int N = vec.size();
    fftw_complex *in, *out;
    fftw_plan p;

    in = fftw_alloc_complex(N);
    out = fftw_alloc_complex(N);

    for (int i = 0; i < N; ++i) {
        in[i][0] = vec[i].real();
        in[i][1] = vec[i].imag();
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    std::vector<std::complex<double>> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = std::complex<double>(out[i][0], out[i][1]);
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return result;
}

std::vector<std::complex<double>> ifft(const std::vector<std::complex<double>>& vec) {
    int N = vec.size();
    fftw_complex *in, *out;
    fftw_plan p;

    in = fftw_alloc_complex(N);
    out = fftw_alloc_complex(N);

    for (int i = 0; i < N; ++i) {
        in[i][0] = vec[i].real();
        in[i][1] = vec[i].imag();
    }

    p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);

    std::vector<std::complex<double>> result(N);
    for (int i = 0; i < N; ++i) {
        result[i] = std::complex<double>(out[i][0], out[i][1]) / static_cast<double>(N);
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    return result;
}

std::vector<std::complex<double>> conjugate(const std::vector<std::complex<double>>& vec) {
    std::vector<std::complex<double>> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = std::conj(vec[i]);
    }
    return result;
}

std::vector<std::complex<double>> elementwise_multiply(const std::vector<std::complex<double>>& vec1, const std::vector<std::complex<double>>& vec2) {
    std::vector<std::complex<double>> result(vec1.size());
    for (size_t i = 0; i < vec1.size(); ++i) {
        result[i] = vec1[i] * vec2[i];
    }
    return result;
}

std::vector<std::complex<double>> flatten(const std::vector<std::vector<std::complex<double>>>& matrix) {
    std::vector<std::complex<double>> result;
    for (const auto& row : matrix) {
        for (const auto& element : row) {
            result.push_back(element);
        }
    }
    return result;
}

std::vector<int> sort_indices(const std::vector<double>& vec, bool ascending = true) {
    std::vector<int> indices(vec.size());
    std::iota(indices.begin(), indices.end(), 0);

    if (ascending) {
        std::sort(indices.begin(), indices.end(), [&vec](int i1, int i2) { return std::abs(vec[i1]) < std::abs(vec[i2]); });
    } else {
        std::sort(indices.begin(), indices.end(), [&vec](int i1, int i2) { return std::abs(vec[i1]) > std::abs(vec[i2]); });
    }

    return indices;
}

std::vector<std::vector<std::complex<double>>> subtract(const std::vector<std::vector<std::complex<double>>>& lhs, const std::vector<std::vector<std::complex<double>>>& rhs) {
    size_t rows = lhs.size();
    size_t cols = lhs[0].size();
    std::vector<std::vector<std::complex<double>>> result(rows, std::vector<std::complex<double>>(cols));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[i][j] = lhs[i][j] - rhs[i][j];
        }
    }

    return result;
}

std::vector<std::complex<double>> subtract_vectors(const std::vector<std::complex<double>> &a, const std::vector<std::complex<double>> &b) {
    assert(a.size() == b.size());

    std::vector<std::complex<double>> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }

    return result;
}


std::vector<double> linspace(double start, double end, int num) {
    std::vector<double> result(num);
    double step = (end - start) / (num - 1);
    for (int i = 0; i < num; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

std::vector<double> subtract_scalar(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] - scalar;
    }
    return result;
}

std::vector<std::vector<std::complex<double>>> multiply_scalar(const std::vector<std::vector<std::complex<double>>>& vec, double scalar) {
    std::vector<std::vector<std::complex<double>>> result(vec.size(), std::vector<std::complex<double>>(vec[0].size()));
    for (size_t i = 0; i < vec.size(); ++i) {
        for (size_t j = 0; j < vec[i].size(); ++j) {
            result[i][j] = vec[i][j] * scalar;
        }
    }
    return result;
}

std::vector<std::complex<double>> multiply_complex_scalar(const std::vector<std::complex<double>>& vec, std::complex<double> scalar) {
    std::vector<std::complex<double>> result(vec.size());

    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }

    return result;
}


std::vector<double> multiply_scalar(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());

    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }

    return result;
}


std::vector<double> abs_vec(const std::vector<std::complex<double>>& vec) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = std::abs(vec[i]);
    }
    return result;
}

double sum(const std::vector<double>& vec) {
    return std::accumulate(vec.begin(), vec.end(), 0.0);
}



std::tuple<std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>> PathParameterEstimator(const std::vector<std::vector<std::complex<double>>>& Y_tr_dd, int DividingNumber, int Ncp, int Nofdm, int NFFT, double fs, double alpha, double beta, double noiseVar, const std::vector<int>& k_subcarrier_p, const std::vector<int>& l_ofdm_p, const std::vector<int>& indx_p_all) {

    // chosen pilot domain receive signal
    std::vector<std::vector<std::complex<double>>> Y_tr_dd_pilot = Y_tr_dd;
    for (int i = 0; i < Y_tr_dd.size(); ++i) {
        for (int j = 0; j < Y_tr_dd[0].size(); ++j) {
            if (std::find(indx_p_all.begin(), indx_p_all.end(), i * Y_tr_dd[0].size() + j) == indx_p_all.end()) {
                Y_tr_dd_pilot[i][j] = std::complex<double>(0, 0);
            }
        }
    }

    int NumDopplerBins = Nofdm;
    int NumDelayBins = NFFT;

    // upsilon
    std::vector<double> kappa = linspace(0, 1 - 1.0 / DividingNumber, DividingNumber);
    std::vector<std::vector<std::complex<double>>> upsilon(NumDopplerBins, std::vector<std::complex<double>>(kappa.size(), std::complex<double>(0, 0)));
    for (int fk = 0; fk < kappa.size(); ++fk) {
        std::vector<double> k = linspace(0, NumDopplerBins - 1, NumDopplerBins);
        std::vector<double> x = subtract_scalar(k, kappa[fk]);
        for (int n = 0; n < NumDopplerBins; ++n) {
            std::vector<std::complex<double>> exp_vec(x.size());
            std::transform(x.begin(), x.end(), exp_vec.begin(), [n, NumDopplerBins](double xi) {
                return std::exp(std::complex<double>(0, 2 * M_PI * (n - 1) * xi / NumDopplerBins));
            });
            for (int i = 0; i < NumDopplerBins; ++i) {
                upsilon[i][fk] += exp_vec[i] / static_cast<double>(NumDopplerBins);
            }
        }
    }

    std::vector<std::vector<std::complex<double>>> Hdd = Y_tr_dd_pilot;
    int M = NFFT;
    int N = Nofdm;
    double SamplingRate = fs;

    std::vector<double> estGains;
    std::vector<double> estDopplers;
    std::vector<double> estDelays;
    std::vector<double> estPhaseOffsets;
    Hdd = multiply_scalar(Hdd, 1 / std::sqrt(M * N));

    kappa = linspace(0, 1 - 1.0 / DividingNumber, DividingNumber);
    int pathindex = 1;
    for (int delay = 0; delay < M; ++delay) {
        std::vector<std::complex<double>> Hprime = Hdd[delay];
        double sumMagnitude = sum(abs_vec(Hprime));
        double estGain = std::numeric_limits<double>::infinity();
        while (true) {

            std::vector<std::vector<std::complex<double>>> crosscorr(kappa.size(), std::vector<std::complex<double>>(N, std::complex<double>(0, 0)));
            for (int idx_kappa = 0; idx_kappa < kappa.size(); ++idx_kappa) {
                crosscorr[idx_kappa] = ifft(elementwise_multiply(fft(Hprime), conjugate(fft(upsilon[idx_kappa]))));
            }

            std::vector<std::complex<double>> crosscorrvec = flatten(crosscorr);

            // Condition 1: find paths until the cross-correlation becomes small
            std::vector<int> largestindices = sort_indices(abs_vec(crosscorrvec), false);
            largestindices.erase(std::remove_if(largestindices.begin(), largestindices.end(), [&](int idx) { return abs(crosscorrvec[idx]) < alpha * sumMagnitude; }), largestindices.end());

            // Condition 2: the cross-correlation is not negligibly small (equally or less than the noise power)
            largestindices.erase(std::remove_if(largestindices.begin(), largestindices.end(), [&](int idx) { return abs(crosscorrvec[idx]) < beta * std::sqrt(noiseVar); }), largestindices.end());

            if (largestindices.empty()) {
                break;
            }

            int largestidx = largestindices[0];

            if (estGain < abs(crosscorrvec[largestidx])) {
                break;
            }
            estGain = abs(crosscorrvec[largestidx]);

            // Estimate Doppler
            double estDoppler = (largestidx - 1) / static_cast<double>(DividingNumber);
            estDoppler = estDoppler - (l_ofdm_p[0] - 1);

            std::complex<double> estDopplerShift = std::complex<double>(cos(2 * M_PI * estDoppler * (Ncp - (delay - (k_subcarrier_p[0] - 1))) / ((M + Ncp) * N)), sin(2 * M_PI * estDoppler * (Ncp - (delay - (k_subcarrier_p[0] - 1))) / ((M + Ncp) * N)));
            std::complex<double> estInitialPhase = crosscorrvec[largestidx] / abs(crosscorrvec[largestidx]) * std::conj(estDopplerShift);

            estGains.push_back(estGain);
            estDopplers.push_back(estDoppler);
            estDelays.push_back(delay - (k_subcarrier_p[0] - 1));
            estPhaseOffsets.push_back(std::arg(estInitialPhase));

            double estDoppler_largest = (largestidx - 1) / static_cast<double>(DividingNumber);
            int intDoppler = static_cast<int>(std::floor(estDoppler_largest));
            int fracDopplerIdx = static_cast<int>(std::floor(std::fmod(estDopplers[pathindex - 1], 1) * DividingNumber + 1 + 1e-9));
            Hprime = subtract_vectors(Hprime, multiply_complex_scalar(upsilon[intDoppler], estGains[pathindex - 1] * std::exp(std::complex<double>(estPhaseOffsets[pathindex - 1]))));

            pathindex++;
        }
    }

    estDopplers = multiply_scalar(estDopplers, 1 / (N * ((M + Ncp) / SamplingRate)));
    estDelays = multiply_scalar(estDelays, 1 / SamplingRate);

    return std::make_tuple(estGains, estDopplers, estDelays, estPhaseOffsets);
}


int main() {
    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1, 1);

    // Generate Y_tr_dd
    std::vector<std::vector<std::complex<double>>> Y_tr_dd(10, std::vector<std::complex<double>>(10));
    for (auto& row : Y_tr_dd) {
        for (auto& element : row) {
            element = std::complex<double>(dis(gen), dis(gen));
        }
    }

    int DividingNumber = 10;
    double fs = 960000;
    double alpha = 0.0200;
    double beta = 0.1000;
    double noiseVar = 1;
    const std::vector<int> k_subcarrier_p = {33};
    const std::vector<int> l_ofdm_p = {8};

    // Generate indx_p_all
    std::vector<int> indx_p_all(10);
    for (int i = 0; i < 10; ++i) {
        indx_p_all[i] = i;
    }

    int NFFT = 64;
    int Nofdm = 14;
    int Ncp = 16;

    auto result = PathParameterEstimator(Y_tr_dd, DividingNumber, Ncp, Nofdm, NFFT, fs, alpha, beta, noiseVar, k_subcarrier_p, l_ofdm_p, indx_p_all);

    std::cout << std::endl;

    return 0;
}
