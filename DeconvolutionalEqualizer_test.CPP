#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using cplx = std::complex<double>;

std::vector<std::vector<cplx>> DeconvolutionalEqualizer(const std::vector<std::vector<cplx>>& Ydd, int DividingNumber, int Ncp, int Nofdm, int NFFT, double fs, const std::string& EqualizationAlgorithm,
    const std::vector<cplx>& pathGains, const std::vector<double>& pathDopplers, const std::vector<double>& pathDelays, const std::vector<double>& pathOffsets, double noiseVar)
{
    int NumDopplerBins = Nofdm;
    int NumDelayBins = NFFT;

    // upsilon
    std::vector<double> kappa(DividingNumber);
    for (int i = 0; i < DividingNumber; ++i) {
        kappa[i] = i / static_cast<double>(DividingNumber);
    }

    std::vector<std::vector<cplx>> upsilon(NumDopplerBins, std::vector<cplx>(DividingNumber, cplx(0, 0)));
    for (int fk = 0; fk < DividingNumber; ++fk) {
        for (int n = 0; n < NumDopplerBins; ++n) {
            std::vector<cplx> temp(NumDopplerBins);
            for (int k = 0; k < NumDopplerBins; ++k) {
                double x = kappa[fk] - k;
                temp[k] = cplx(cos(2 * M_PI * (n) * x / NumDopplerBins), sin(2 * M_PI * (n) * x / NumDopplerBins));
            }
            cplx sum = std::accumulate(temp.begin(), temp.end(), cplx(0, 0));
            upsilon[n][fk] = sum / static_cast<double>(NumDopplerBins);
        }
    }

    int NumSubcarriers = NFFT;
    int NumSymbols = Nofdm;

    int P = pathGains.size();
    int M = NumSubcarriers;
    int N = NumSymbols;
    double SamplingRate = fs;

    // Convert delay and Doppler into grid indices
    std::vector<double> pathDelayInSample(P);
    std::vector<double> pathDopplerInSample(P);
    for (int p = 0; p < P; ++p) {
        pathDelayInSample[p] = pathDelays[p] * SamplingRate;
        pathDopplerInSample[p] = pathDopplers[p] * N * ((M + Ncp) / SamplingRate);
    }

    /// Calc Lambda (an effective channel matrix)
    std::vector<std::vector<std::vector<cplx>>> lambda(P, std::vector<std::vector<cplx>>(M, std::vector<cplx>(N, cplx(0, 0))));
    for (int p = 0; p < P; ++p) {
        for (int l = 0; l < M; ++l) {
            int intDoppler = static_cast<int>(round(pathDopplerInSample[p])) % N;
            int fracDopplerIndex = static_cast<int>(round((pathDopplerInSample[p] - intDoppler) * DividingNumber)) % DividingNumber;

            cplx hp = pathGains[p] * (pathDelayInSample[p] == l);
            for (int n = 0; n < N; ++n) {
                lambda[p][l][n] = hp * upsilon[intDoppler][fracDopplerIndex];
                intDoppler = (intDoppler + 1) % NumDopplerBins;
            }
        }
    }

    std::vector<std::vector<cplx>> Ydd_dft = fft2(Ydd);
    std::vector<std::vector<cplx>> Xdd(M, std::vector<cplx>(N, cplx(0, 0)));

    for (int l = 0; l < M; ++l) {
        std::vector<cplx> psi(P);
        std::vector<cplx> phi(P);
        for (int p = 0; p < P; ++p) {
            psi[p] = cplx(cos(2 * M_PI * pathDopplerInSample[p] * (Ncp - pathDelayInSample[p] + l) / ((M + Ncp) * N)), sin(2 * M_PI * pathDopplerInSample[p] * (Ncp - pathDelayInSample[p] + l) / ((M + Ncp) * N)));
            phi[p] = cplx(cos(pathOffsets[p]), sin(pathOffsets[p]));
        }

        std::vector<std::vector<cplx>> Hdd(M, std::vector<cplx>(N, cplx(0, 0)));
        for (int m = 0; m < M; ++m) {
            for (int n = 0; n < N; ++n) {
                for (int p = 0; p < P; ++p) {
                    Hdd[m][n] += psi[p] * phi[p] * lambda[p][m][n];
                }
            }
        }

        std::vector<std::vector<cplx>> Hdd_dft = fft2(Hdd);

        std::vector<std::vector<cplx>> Xdd_tmp;
        if (EqualizationAlgorithm == "ZF") {
            Xdd_tmp = ifft2(elementwise_divide(Ydd_dft, Hdd_dft));
        } else if (EqualizationAlgorithm == "MMSE") {
            std::vector<std::vector<cplx>> eq = elementwise_divide(conjugate(Hdd_dft), elementwise_add(elementwise_multiply(Hdd_dft, conjugate(Hdd_dft)), noiseVar));
            Xdd_tmp = ifft2(elementwise_multiply(Ydd_dft, eq));
        } else {
            throw std::runtime_error("Invalid equalization algorithm specified.");
        }

        for (int n = 0; n < N; ++n) {
            Xdd[l][n] = Xdd_tmp[l][n];
        }
    }

    return Xdd;
}