#include <iostream>
#include <complex>
#include <vector>
#include <cmath>
#include <fstream>


using namespace std;

vector<complex<double>> add_CFO_myself(const vector<complex<double>>& y, double freq_offset, double fs) {
    vector<complex<double>> y_CFO(y.size());

    for(int i=0; i<y.size(); i++) {
        double time_base = static_cast<double>(i)/fs;
        complex<double> phase_rotation = exp(complex<double>(0, 1)*2.0*M_PI*static_cast<double>(freq_offset)*time_base);
        y_CFO[i] = y[i] * phase_rotation;
    }

    return y_CFO;
}

int main() {

    // Open input file
    ifstream inputFile("y_1120.txt");
    if(!inputFile.is_open()) {
        cerr << "Failed to open input file." << endl;
        return 1;
    }

    // Read data from input file
    vector<double> y_real;
    double value;
    while(inputFile >> value) {
        y_real.push_back(value);
    }

    // Convert real data to complex numbers
    vector<complex<double>> y(y_real.size());
    for(int i=0; i<y_real.size(); i++) {
        y[i] = complex<double>(y_real[i], 0);
    }

    // Close input file
    inputFile.close();

    
    cout << "import the number of y_real :" << y_real.size() << endl;
    
    double freq_offset = 2.571428571428572e+03;
    double fs = 960000;
    
    vector<complex<double>> y_CFO=add_CFO_myself(y,freq_offset, fs);


    for (int n = 0; n < y_CFO.size(); ++n) {
        std::cout << y_CFO[n] << " ";
    }
    std::cout << std::endl;
    std::cout << "output, the number of y_CFO , Output size y_CFO : " << y_CFO.size() << std::endl;

    return 0;
}
