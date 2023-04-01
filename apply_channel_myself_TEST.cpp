#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

vector<complex<double>> apply_channel_myself(const vector<double>& x, double fs, const vector<vector<complex<double>>>& chan_coef, const vector<int>& delay_taps, const vector<double>& Dopplers, int taps, double startTime) {
    
    int Lx = x.size();
    int Lg = *max_element(delay_taps.begin(), delay_taps.end());
    int L = Lx + Lg;

    vector<vector<complex<double>>> y(taps, vector<complex<double>>(L));
    for (int p = 0; p < taps; ++p) {
        int l = delay_taps[p] + 1;
        for (int n = 0; n < Lx; ++n) {
            y[p][l - 1 + n] = x[n] * (chan_coef[p][0] * exp(complex<double>(0, 1) * 2.0 * M_PI / fs * Dopplers[p] * (n + startTime)));
        }
    }

    vector<complex<double>> z(L);
    for (int n = 0; n < L; ++n) {
        complex<double> sum(0, 0);
        for (int p = 0; p < taps; ++p) {
            sum += y[p][n];
        }
        z[n] = sum;
    }

    return z;
}

int main() {
     /*std::vector<std::vector<double>> x(1120); 
    
    
     for (int i = 0; i < 1120; i++) {
        std::ifstream infile("vector.txt");
        std::vector<double> temp;
        double value;
        while (infile >> value) {
            temp.push_back(value);
        }
        
        
        x[i] = temp;

    }
    
        std::vector<double> flattened;
        for (const auto& row : x) {
        flattened.insert(flattened.end(), row.begin(), row.end());
    }*/
    

    
    /*ifstream file("S_TR1120.txt");
    vector<double> x;

    double num;
    while (file >> num) {
    x.push_back(num);
    }

    for (double val : x) {
        cout << val << endl;
    }*/

    // 打开文件
    ifstream inputFile("vector.txt");

    // 如果文件无法打开，打印错误消息并退出程序
    if (!inputFile)
    {
        cerr << "File kaputt" << endl;
        return 1;
    }

    // 读取文件中的数据
    vector<double> x;
    double value;
    while (inputFile >> value)
    {
        x.push_back(value);
    }

    // 关闭文件
    inputFile.close();

    // 后面的代码可以使用变量 x，例如：
    cout << "import the number of x :" << x.size() << endl;


    double fs = 960000.0;
    std::vector<std::vector<std::complex<double>>> chan_coef = {
        {0.186718629470780 + 0.961757653745168i},{0.535853729921359 - 0.394431464403843i},{-0.667668589267060 + 0.897060842982048i},{0.197819334266011 + 0.166438670117188i},{0.103310654628862 - 0.020435855018713i},{0.087061071632880i},{-0.067259779514184 - 0.031794800595789i},{0.029887702296789 - 0.010829318054259i},{0.177567168755591 + 0.073921756298147i}
    };
    
    std::vector<int> delay_taps = {0,0,0,0,0,1,1,2,2};
    std::vector<double> Dopplers = {97.140861087160470, 3.586870755326057e+02, -2.068942837778906e+02, 3.613357363967352e+02, 2.162053363231923e+02, 3.392065461847558e+02, -1.604653304907923e+02, 18.017467012018810, -15.987573761070944};
    int taps = 9;
    double startTime = 0.0;

    std::vector<std::complex<double>> z = apply_channel_myself(x, fs, chan_coef, delay_taps, Dopplers, taps, startTime);

    /*for (size_t n = 0; n < z.size(); ++n) {
        cout << z[n].real() << " + " << z[n].imag() << "i ";
    }
    cout << endl;

    cout << "Output size x : " << x.size() << endl;

    // 获取向量的大小
    std::vector<complex<double>>::size_type size = z.size();
    std::cout << "Output size z :" << size << std::endl;*/

    for (int n = 0; n < z.size(); ++n) {
        std::cout << z[n] << " ";
    }
    std::cout << std::endl;
    std::cout << "output, the number of z , Output size z : " << z.size() << std::endl;

    return 0;

    /*// 测试代码
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
    double fs = 44100.0;
    std::vector<std::vector<std::complex<double>>> chan_coef = {{1.0 + 1.0i, 2.0 + 2.0i, 3.0 + 3.0i},
                                                                {4.0 + 4.0i, 5.0 + 5.0i, 6.0 + 6.0i}};
    std::vector<int> delay_taps = {0, 1};
    std::vector<double> Dopplers = {1000.0, 2000.0};
    int taps = 2;
    double startTime = 0.0;

    std::vector<std::complex<double>> z = apply_channel_myself(x, fs, chan_coef, delay_taps, Dopplers, taps, startTime);

    for (int n = 0; n < z.size(); ++n) {
        std::cout << z[n] << " ";
    }
    std::cout << std::endl;
    std::cout << "output, the number of z , Output size z : " << z.size() << std::endl;

    return 0;*/

}


