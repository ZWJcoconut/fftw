#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <complex>
#include <Eigen/Dense>
#include <fftw3.h>

using namespace std;
using namespace Eigen;

typedef complex<double> cd;

void sfft(MatrixXcd &Y_tr_ft) {
    // Implement sFFT function here
}

double CFO_estimation_afterFFT(const VectorXcd &rxSig_tr,
                               const vector<int> &taps_subcarrier_location,
                               int k_subcarrier_p,
                               int l_ofdm_p,
                               double fs,
                               int Ncp,
                               double noiseVar,
                               int NFFT,
                               int Nofdm,
                               const vector<int> &indx_p_all) {
    double est_CFO = 0;

    // S/P
    MatrixXcd rxSig_tr_shaped(NFFT + Ncp, Nofdm);
    for (int i = 0; i < Nofdm; i++) {
        for (int j = 0; j < NFFT + Ncp; j++) {
            rxSig_tr_shaped(j, i) = rxSig_tr(i * (NFFT + Ncp) + j);
        }
    }

    // Remove CP
    rxSig_tr_shaped = rxSig_tr_shaped.bottomRows(NFFT);

    // FFT or Wigner transform
    MatrixXcd Y_tr_ft(NFFT, Nofdm);
    fftw_complex *in, *out;
    fftw_plan p;
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NFFT);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * NFFT);
    p = fftw_plan_dft_1d(NFFT, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (int i = 0; i < Nofdm; i++) {
        for (int j = 0; j < NFFT; j++) {
            in[j][0] = rxSig_tr_shaped(j, i).real();
            in[j][1] = rxSig_tr_shaped(j, i).imag();
        }

        fftw_execute(p);

        for (int j = 0; j < NFFT; j++) {
            Y_tr_ft(j, i) = cd(out[j][0], out[j][1]) / sqrt(NFFT);
        }
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);

    // sFFT
    sfft(Y_tr_ft);

    // Pilot grid
    MatrixXcd Y_tr_dd_pilot = MatrixXcd::Zero(NFFT, Nofdm);
    for (const auto &idx : indx_p_all) {
        int row = idx % NFFT;
        int col = idx / NFFT;
        Y_tr_dd_pilot(row, col) = Y_tr_ft(row, col);
    }

    MatrixXcd Y_tr_dd = Y_tr_dd_pilot;

    int NumDopplerBins = Nofdm;
    int NumDelayBins = NFFT;
    int DividingNumber = 10;

    // upsilon
    vector<double> kappa;
    for (int i = 0; i < DividingNumber; i++) {
        kappa.push_back(i / static_cast<double>(DividingNumber));
    }

    MatrixXcd upsilon(NumDopplerBins, kappa.size());
    for (size_t fk = 0; fk < kappa.size(); fk++) {
        for (int n = 1; n <= NumDopplerBins; n++) {
            VectorXcd x = kappa[fk] * VectorXd::LinSpaced(NumDopplerBins, 0, NumDopplerBins - 1);
            upsilon.col(fk) += (1.0 / NumDopplerBins) * (VectorXcd::Unit(NumDopplerBins, n - 1) * x.array().exp().transpose());
        }
    }

    MatrixXcd Hdd = Y_tr_dd;
    int M = NFFT;
    int N = Nofdm;
    double SamplingRate = fs;

    // Parameters for choosing largest taps
    double alpha = 1.0 / 50;
    double beta = 1.0 / 10;

    // Variable initialization
    vector<double> estGains;
    vector<double> estDopplers;
    vector<double> estDelays;
    vector<double> estPhaseOffsets;
    Hdd = Hdd / sqrt(M * N);

    int pathindex = 1;
    for (int p = 0; p < taps_subcarrier_location.size(); p++) {
        int delay = taps_subcarrier_location[p] - 1;

        VectorXcd Hprime = Hdd.row(delay);
        double sumMagnitude = Hprime.cwiseAbs().sum();
        double estGain = INFINITY;

        while (true) {
            MatrixXcd crosscorr(kappa.size(), N);
            for (size_t idx_kappa = 0; idx_kappa < kappa.size(); idx_kappa++) {
                crosscorr.row(idx_kappa) = (Hprime.array() * (upsilon.col(idx_kappa).array().conjugate())).matrix().fft();
            }

            VectorXcd crosscorrvec = Map<VectorXcd>(crosscorr.data(), crosscorr.size());
            vector<pair<double, int>> largest_indices;
            for (int i = 0; i < crosscorrvec.size(); i++) {
                largest_indices.push_back(make_pair(abs(crosscorrvec(i)), i));
            }

            sort(largest_indices.begin(), largest_indices.end(), greater<pair<double, int>>());

            while (!largest_indices.empty() && largest_indices.back().first < alpha * sumMagnitude) {
                largest_indices.pop_back();
            }

            while (!largest_indices.empty() && abs(crosscorrvec(largest_indices.back().second)) < beta * sqrt(noiseVar)) {
                largest_indices.pop_back();
            }

            if (largest_indices.empty()) {
                break;
            }

            int largestidx = largest_indices[0].second;

            if (estGain < abs(crosscorrvec(largestidx))) {
                break;
            }
            estGain = abs(crosscorrvec(largestidx));

            double estDoppler = (largestidx - 1) / static_cast<double>(DividingNumber);
            estDoppler = estDoppler - (l_ofdm_p - 1);

            cd estDopplerShift = exp(cd(0, 1) * 2.0 * M_PI * estDoppler * (Ncp - (delay - (k_subcarrier_p - 1))) / ((M + Ncp) * N));
            cd estInitialPhase = crosscorrvec(largestidx) / abs(crosscorrvec(largestidx)) * pow(estDopplerShift, -1);
            
            estGains.push_back(estGain);
            estDopplers.push_back(estDoppler);
            estDelays.push_back(delay - (k_subcarrier_p - 1));
            estPhaseOffsets.push_back(arg(estInitialPhase));

            double estDoppler_largest = (largestidx - 1) / static_cast<double>(DividingNumber);
            int intDoppler = floor(estDoppler_largest);
            int fracDopplerIdx = floor(mod(estDopplers[pathindex - 1], 1) * DividingNumber + 1 + 1e-9);

            Hprime -= estGains[pathindex - 1] * exp(cd(0, estPhaseOffsets[pathindex - 1])) * circshift(upsilon.col(fracDopplerIdx), intDoppler).transpose();
            pathindex++;
        }
    }

    vector<double> estDelays_unique = unique(estDelays);
    vector<double> estDopplers_largest(estDelays_unique.size(), 0);
    for (size_t i = 0; i < estDelays_unique.size(); i++) {
        vector<double> estGains_idx, estDopplers_idx;
        for (size_t j = 0; j < estDelays.size(); j++) {
            if (estDelays[j] == estDelays_unique[i]) {
                estGains_idx.push_back(estGains[j]);
                estDopplers_idx.push_back(estDopplers[j]);
            }
        }
        auto max_gain = max_element(estGains_idx.begin(), estGains_idx.end());
        size_t idx_largest_estGain = distance(estGains_idx.begin(), max_gain);
        estDopplers_largest[i] = estDopplers_idx[idx_largest_estGain];
    }

    double est_CFO = estDopplers_largest[0];
    return est_CFO;
}
