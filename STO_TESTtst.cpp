#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <fstream>

std::pair<int, std::vector<int>> STO_estimation_afterFFT(const std::vector<std::complex<double>> &rxSig_tr, int NFFT, int Nofdm, int k_subcarrier_p, int l_ofdm_p, const std::vector<int> &k_subcarrier_p_all, const std::vector<int> &l_ofdm_p_all, int Ncp) {

    int rows = NFFT + Ncp;
    std::vector<std::vector<std::complex<double>>> rxSig_tr_shaped(rows, std::vector<std::complex<double>>(Nofdm));

    // Reshape rxSig_tr
    for (int i = 0; i < Nofdm; i++) {
        for (int j = 0; j < rows; j++) {
            rxSig_tr_shaped[j][i] = rxSig_tr[i * rows + j];
        }
    }

    // Remove CP
    rxSig_tr_shaped.erase(rxSig_tr_shaped.begin(), rxSig_tr_shaped.begin() + Ncp);

    std::vector<std::vector<std::complex<double>>> yp_mat(k_subcarrier_p_all.size(), std::vector<std::complex<double>>(l_ofdm_p_all.size()));

    for (int i = 0; i < k_subcarrier_p_all.size(); i++) {
        for (int j = 0; j < l_ofdm_p_all.size(); j++) {
            yp_mat[i][j] = rxSig_tr_shaped[k_subcarrier_p_all[i] - 1][l_ofdm_p_all[j] - 1];
        }
    }

    int sliding_win = 1;
    std::vector<std::vector<std::complex<double>>> P_mat(k_subcarrier_p_all.size(), std::vector<std::complex<double>>(l_ofdm_p_all.size() - sliding_win));

    for (int idx_k = 0; idx_k < k_subcarrier_p_all.size(); idx_k++) {
        std::vector<std::complex<double>> r_vec = yp_mat[idx_k];
        for (int l = 0; l < l_ofdm_p_all.size() - sliding_win; l++) {
            std::vector<std::complex<double>> r_part1(r_vec.begin() + l, r_vec.begin() + l + sliding_win);
            std::vector<std::complex<double>> r_part2(r_vec.begin() + l + 1, r_vec.begin() + l + 1 + sliding_win);

            std::complex<double> P_sum(0.0, 0.0);
            for (int i = 0; i < r_part1.size(); i++) {
                P_sum += std::conj(r_part1[i]) * r_part2[i];
            }
            P_mat[idx_k][l] = P_sum;
        }
    }

    std::vector<std::complex<double>> Pd(k_subcarrier_p_all.size());
    for (int i = 0; i < k_subcarrier_p_all.size(); i++) {
        Pd[i] = std::accumulate(P_mat[i].begin(), P_mat[i].end(), std::complex<double>(0.0, 0.0));
    }

    std::vector<double> Pd_abs(k_subcarrier_p_all.size());
    for (int i = 0; i < k_subcarrier_p_all.size(); i++) {
        Pd_abs[i] = std::abs(Pd[i]);
    }

    std::vector<double> Pd_abs_sort = Pd_abs;
    std::vector<size_t> largestindices(Pd_abs.size());
    std::iota(largestindices.begin(), largestindices.end(), 0);

    std::sort(largestindices.begin(), largestindices.end(), [&Pd_abs_sort](size_t i1, size_t i2) {
        return Pd_abs_sort[i1] > Pd_abs_sort[i2];
    });

    std::sort(Pd_abs_sort.rbegin(), Pd_abs_sort.rend());

    int max_path = 10;
    int Num_noise = Pd_abs.size() - max_path;
    std::vector<double> Pd_abs_noise(Pd_abs_sort.end() - Num_noise, Pd_abs_sort.end());
    double Pd_abs_noise_max = *std::max_element(Pd_abs_noise.begin(), Pd_abs_noise.end());

    // limit1
    double limit1 = 4 * Pd_abs_noise_max;

    // limit2
    double sumMagnitude = std::abs(std::accumulate(Pd.begin(), Pd.end(), std::complex<double>(0.0, 0.0)));
    double threshAlpha = 1.0 / 50.0;
    double limit2 = threshAlpha * sumMagnitude;

    double limit = std::max(limit1, limit2);

    std::vector<int> index_taps_all;
    for (size_t i = 0; i < Pd_abs.size(); i++) {
        if (std::abs(Pd[i]) > limit) {
            index_taps_all.push_back(i);
        }
    }

    int idx_k_first_tap = index_taps_all[0];
    int first_taps_k = k_subcarrier_p_all[idx_k_first_tap];

    int theta_d = k_subcarrier_p - first_taps_k;

    int est_nSTO = theta_d;

    std::vector<int> taps_subcarrier_location(index_taps_all.size());
    for (size_t i = 0; i < index_taps_all.size(); i++) {
        taps_subcarrier_location[i] = k_subcarrier_p_all[index_taps_all[i]];
    }

    std::vector<int> first_taps_subcarrier_location = taps_subcarrier_location;

    return std::make_pair(est_nSTO, first_taps_subcarrier_location);
}


int main() {
    // open a file
    std::ifstream inputFile("rxSig_tr_1120.txt");

    // If the file cannot be opened, print an error message and exit the program
    if (!inputFile)
    {
        std::cerr << "File kaputt" << std::endl;
        return 1;
    }

    // read data from file
    std::vector<std::complex<double>> rxSig_tr;
    double value;
    while (inputFile >> value)
    {
        double real = value;
        inputFile >> value;
        double imag = value;
        rxSig_tr.push_back(std::complex<double>(real, imag));
    }

    // close the file
    inputFile.close();

    int NFFT = 64;
    int Nofdm = 14;
    int k_subcarrier_p = 33;
    int l_ofdm_p = 8;
    std::vector<int> k_subcarrier_p_all = {8,9,10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58};
    std::vector<int> l_ofdm_p_all = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
    int Ncp = 16;

    auto result = STO_estimation_afterFFT(rxSig_tr, NFFT, Nofdm, k_subcarrier_p, l_ofdm_p, k_subcarrier_p_all, l_ofdm_p_all, Ncp);
    int est_nSTO = result.first;
    std::vector<int> first_taps_subcarrier_location = result.second;

    std::cout << "est_nSTO: " << est_nSTO << std::endl;
    std::cout << "first_taps_subcarrier_location: ";
    for (auto loc : first_taps_subcarrier_location) {
        std::cout << loc << " ";
    }
    std::cout << std::endl;

    return 0;
}