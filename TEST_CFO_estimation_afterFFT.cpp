#include <iostream>
#include <complex>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

const double pi = 3.14159265358979323846;

void circshift(std::vector<int>& v, int k) {
    int n = v.size();
    k = (k % n + n) % n; // 处理 k 取负数的情况
    std::rotate(v.rbegin(), v.rbegin() + k, v.rend());
}


vector<double> linspace(double start, double end, int numPoints)
{
    vector<double> arr(numPoints);
    double step = (end - start) / (numPoints - 1);
    for (int i = 0; i < numPoints; i++)
        arr[i] = start + i * step;
    return arr;
}

vector<vector<complex<double>>> sfft(vector<vector<complex<double>>> x)
{
    int M = x.size();
    int N = x[0].size();
    vector<vector<complex<double>>> y(M, vector<complex<double>>(N));

    for (int k = 0; k < N; k++)
    {
        vector<complex<double>> x_k(M);
        for (int m = 0; m < M; m++)
            x_k[m] = x[m][k];

        vector<complex<double>> y_k(M);
        for (int n = 0; n < M; n++)
        {
            vector<double> indices = linspace(0, M - 1, M);
            vector<complex<double>> x_k_shifted(M);
            for (int m = 0; m < M; m++)
            {
                int ind = fmod(m - n + M, M);
                x_k_shifted[m] = x_k[ind];
            }
            std::complex<double> sum = 0;
            for (int m = 0; m < M; m++){
                sum += x_k_shifted[m] * exp(-1i * 2.0 * M_PI * m * n / M);
                y_k[n] = sum;
            }
        }

        for (int m = 0; m < M; m++)
            y[m][k] = y_k[m];
    }



    return y;
}

vector<double> CFO_estimation_afterFFT(vector<double> rxSig_tr, vector<int> taps_subcarrier_location, int k_subcarrier_p, int l_ofdm_p, double fs, int Ncp, double noiseVar, int NFFT, int Nofdm, vector<vector<int>> indx_p_all)
{
    vector<double> est_CFO;

    vector<vector<complex<double>>> rxSig_tr_shaped(NFFT + Ncp, vector<complex<double>>(Nofdm));
    for (int i = 0; i < NFFT + Ncp; i++)
        for (int j = 0; j < Nofdm; j++)
            rxSig_tr_shaped[i][j] = rxSig_tr[i * Nofdm + j];

    // Remove CP
    rxSig_tr_shaped.erase(rxSig_tr_shaped.begin(), rxSig_tr_shaped.begin() + Ncp);

    // Function to perform the FFT on the input signal
    vector<complex<double>> fft(vector<complex<double>>& input_signal) 
        int N = input_signal.size();
        vector<complex<double>> output_signal(N);

        for (int k = 0; k < N; k++) {
            complex<double> sum(0.0, 0.0);
            for (int n = 0; n < N; n++) {
                double angle = 2.0 * M_PI * k * n / N;
                sum += input_signal[n] * complex<double>(cos(angle), -sin(angle));
            }
            output_signal[k] = sum;
        }

    

    // Function to perform the SFFT on the input signal
    vector<complex<double>> sfft(vector<complex<double>>& input_signal) 
        int N = input_signal.size();
        vector<complex<double>> output_signal(N);

        for (int k = 0; k < N; k++) {
            complex<double> sum(0.0, 0.0);
            for (int m = 0; m < N; m++) {
                int n = m + k;
                if (n >= N) {
                    n -= N;
                }
                sum += input_signal[m] * conj(input_signal[n]);
            }
            output_signal[k] = sum / sqrt(N);
        }

    



    vector<complex<double>> rxSig_tr_shaped(NFFT);

    // Perform FFT on input signal
    vector<complex<double>> Y_tr_ft = fft(rxSig_tr_shaped);

    // Perform SFFT on FFT output signal
    vector<complex<double>> Y_tr_dd = sfft(Y_tr_ft);


    // pilot grid, take out received signal in pilot domain and set zero outside pilot domain
    vector<vector<complex<double>>> Y_tr_dd_pilot(Y_tr_dd.size(), vector<complex<double>>(Y_tr_dd[0].size(), complex<double>(0, 0)));
    for (int i = 0; i < indx_p_all.size(); i++) {
        int row = indx_p_all[i] / Y_tr_dd[0].size();
        int col = indx_p_all[i] % Y_tr_dd[0].size();
        Y_tr_dd_pilot[row][col] = Y_tr_dd[row][col];
    }
        Y_tr_dd = Y_tr_dd_pilot;

        int NFFT = Y_tr_dd.size();
        int Nofdm = Y_tr_dd[0].size();
        int NumDopplerBins = Nofdm;
        int NumDelayBins = NFFT;
        int DividingNumber = 10;

    // upsilon
    vector<double> kappa(DividingNumber);   // Define kappa
    for (int i = 0; i < DividingNumber; i++) {
        kappa[i] = (double)i / DividingNumber;
    }
    vector<vector<complex<double>>> upsilon(NumDopplerBins, vector<complex<double>>(kappa.size(), complex<double>(0, 0)));  // Define upsilon
    for (int fk = 0; fk < kappa.size(); fk++) {
            vector<int> k(NumDopplerBins);
        for (int i = 0; i < NumDopplerBins; i++) {
            k[i] = i;
        }

            for (int n = 0; n < NumDopplerBins; n++) {
                complex<double> sum = 0;
                for (int i = 0; i < NumDopplerBins; i++) {
                    double x = kappa[fk] - k[i];
                    sum += 1.0 / NumDopplerBins * std::exp(std::complex<double> (0, 1) * 2 * M_PI * (n * i) / NumDopplerBins * x);
                }
                    upsilon[n][fk] = sum;
            }
    }

    vector<vector<complex<double>>> Hdd = Y_tr_dd;
    int M = NFFT;
    int N = Nofdm;
    double SamplingRate = fs;

    // parameter for chosing largest taps
    double alpha = 1.0 / 50;
    double beta = 1.0 / 10;

    // variable initalization
    vector<double> estGains;
    vector<double> estDopplers;
    vector<double> estDelays;
    vector<double> estPhaseOffsets;
    for (int i = 0; i < Hdd.size(); i++) {
        for (int j = 0; j < Hdd[0].size(); j++) {
            Hdd[i][j] /= sqrt(M * N);
        }
    }

    //fractional dopplers between integer dopplers

    int pathindex = 0;
    for(int p=0; p<taps_subcarrier_location.size(); ++p)
    {
        int delay = taps_subcarrier_location[p]-1;

        std::vector<std::complex<double>> Hprime(Nofdm);
        for(int i=0; i<Nofdm; ++i)
        {
            Hprime[i] = Hdd[delay+1][i];
        }

        double sumMagnitude = 0.0;
        for(int i=0; i<Nofdm; ++i)
        {
            sumMagnitude += std::abs(Hprime[i]);
        }

        double estGain = std::numeric_limits<double>::infinity();
        while(true)
        {
            // Take cross-correlation b/w Hdd and Upsilon (see above upsilon expression)
            std::vector<std::vector<std::complex<double>>> crosscorr(DividingNumber, std::vector<std::complex<double>>(Nofdm));
            for(int idx_kappa=0; idx_kappa<DividingNumber; ++idx_kappa)
            {
                for(int n=0; n<Nofdm; ++n)
                {
                    std::complex<double> sum = 0.0;
                    for(int k=0; k<Nofdm; ++k)
                    {
                        sum += std::exp(std::complex<double>(0, -2*M_PI*k*(kappa[idx_kappa]-n)/Nofdm)) * Hprime[k];
                    }
                    crosscorr[idx_kappa][n] = sum/static_cast<double>(Nofdm);
                }
            }

            std::vector<std::complex<double>> crosscorrvec(DividingNumber*Nofdm);
            for(int i=0; i<DividingNumber; ++i)
            {
                for(int j=0; j<Nofdm; ++j)
                {
                    crosscorrvec[i*Nofdm+j] = crosscorr[i][j];
                }
            }

            // Condition 1: find paths until the cross-correlation becomes small
            std::vector<int> largestindices;
            for(int i=0; i<crosscorrvec.size(); ++i)
            {
                if(std::abs(crosscorrvec[i]) >= alpha*sumMagnitude)
                {
                    largestindices.push_back(i);
                }
            }

            // Condition 2: the cross-correlation is not negligibly small (equally or less than the noise power)
            for(auto it=largestindices.begin(); it!=largestindices.end();)
            {
                if(std::abs(crosscorrvec[*it]) < beta*std::sqrt(noiseVar))
                {
                    it = largestindices.erase(it);
                }
                else
                {
                    ++it;
                }
            }

            if(largestindices.empty())
            {
                break;
            }

            int largestidx = largestindices[0];
            if(estGain < std::abs(crosscorrvec[largestidx]))
            {
                break;
            }
            double estGain = abs(crosscorrvec[largestidx]);

            double estDoppler = (largestidx-1) / DividingNumber;
            estDoppler = estDoppler - (l_ofdm_p-1);

            double estDopplerShift = exp(1i * 2 * M_PI * estDoppler * (Ncp - (delay - (k_subcarrier_p-1))) / ((M+Ncp)*N));
            std::complex<double> estInitialPhase = crosscorrvec[largestidx] / std::abs(crosscorrvec[largestidx]) * pow(estDopplerShift, -1);

            estGains.push_back(estGain);
            estDopplers.push_back(estDoppler);
            estDelays.push_back(delay - (k_subcarrier_p-1));
            estPhaseOffsets.push_back(arg(estInitialPhase));

            double estDoppler_largest = (largestidx-1) / DividingNumber;
            int intDoppler = floor(estDoppler_largest);
            int fracDopplerIdx = floor(fmod(estDopplers[pathindex], 1) * DividingNumber + 1 + 1e-9);

            // remove current tap from received current channel response(Hprime) for next tap
            // Hprime = Hprime - estGains[pathindex] * exp(1i * estPhaseOffsets[pathindex]) * circshift(upsilon.col(fracDopplerIdx), intDoppler).adjoint();
            Hprime = Hprime - estGains[pathindex] * std::exp(std::complex<double>(0, 1) * estPhaseOffsets[pathindex]) * circshift(upsilon, -intDoppler)[fracDopplerIdx].adjoint();


            pathindex++;
        }
    }

    vector<int> estDelays_unique;
    for(int i = 0; i < estDelays.size(); i++){
        if(find(estDelays_unique.begin(), estDelays_unique.end(), estDelays[i]) == estDelays_unique.end()){
            estDelays_unique.push_back(estDelays[i]);
        }
    }

    vector<double> estDopplers_largest(estDelays_unique.size(), 0.0);
    for(int i = 0; i < estDelays_unique.size(); i++){
        vector<int> idx;
        for(int j = 0; j < estDelays.size(); j++){
            if(estDelays[j] == estDelays_unique[i]){
                idx.push_back(j);
            }
        }
        vector<double> estGains_idx;
        vector<double> estDopplers_idx;
        for(int j = 0; j < idx.size(); j++){
            estGains_idx.push_back(estGains[idx[j]]);
            estDopplers_idx.push_back(estDopplers[idx[j]]);
        }
        vector<int> idx_largest_estGain;
        double max_estGains_idx = *max_element(estGains_idx.begin(), estGains_idx.end());
        for(int j = 0; j < estGains_idx.size(); j++){
            if(estGains_idx[j] == max_estGains_idx){
                idx_largest_estGain.push_back(j);
            }
        }
        estDopplers_largest[i] = estDopplers_idx[idx_largest_estGain[0]];
    }

    double est_CFO = estDopplers_largest[0];


    return est_CFO;
}


int main(){
    
    // 打开文件
    ifstream inputFile("rxSig_tr_1120.txt");

    // 如果文件无法打开，打印错误消息并退出程序
    if (!inputFile)
    {
        cerr << "File kaputt" << endl;
        return 1;
    }

    // 读取文件中的数据
    vector<complex<double>> rxSig_tr;
    double value;
    while (inputFile >> value)
    {
        rxSig_tr.push_back(value);
    }

    // 关闭文件
    inputFile.close();

    // 后面的代码可以使用变量 x，例如：
    cout << "import the number of rxSig_tr :" << rxSig_tr.size() << endl;

    
    vector<int> taps_subcarrier_location = {33}; 
    int k_subcarrier_p = 33; 
    int l_ofdm_p = 8; 
    double fs = 960000.0; 
    double Ncp = 16; 
    double noiseVar = 1; 
    int NFFT = 64; 
    int Nofdm = 14;

    // 打开文件
    ifstream inputFile_indx_p_all_663("indx_p_all_663.txt");

    // 如果文件无法打开，打印错误消息并退出程序
    if (!inputFile_indx_p_all_663)
    {
        cerr << "File kaputt" << endl;
        return 1;
    }

    // 读取文件中的数据
    vector<int> indx_p_all;
    double value_indx_p_all;
    while (inputFile_indx_p_all_663 >> value_indx_p_all)
    {
        indx_p_all.push_back(value_indx_p_all);
    }

    // 关闭文件
    inputFile.close();

    // 后面的代码可以使用变量 x，例如：
    cout << "import the number of indx_p_all :" << indx_p_all.size() << endl;

    

    vector<double> est_CFO = CFO_estimation_afterFFT(rxSig_tr,taps_subcarrier_location,k_subcarrier_p,l_ofdm_p,fs,Ncp,noiseVar,NFFT,Nofdm,indx_p_all)



    cout << "import the number of rxSig_tr :" << rxSig_tr.size() << endl;
    cout << "import the number of indx_p_all :" << indx_p_all.size() << endl;
    std::cout << std::endl;
    std::cout << "output, the number of z , Output size z : " << est_CFO.size() << std::endl;

    return 0;
}


