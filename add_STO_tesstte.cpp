#include <complex>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
using namespace std;

std::vector<std::complex<double>> add_STO(std::vector<std::complex<double>> y, int iSTO) {
    std::vector<std::complex<double>> y_STO;
    int N = y.size();
    
    if (iSTO >= 0) {
        y_STO.resize(N);
        for (int i = 0; i < N - iSTO; i++) {
            y_STO[i] = y[i + iSTO];
        }
        for (int i = N - iSTO; i < N; i++) {
            y_STO[i] = 0;
        }
    } else {
        y_STO.resize(N);
        for (int i = 0; i < -iSTO; i++) {
            y_STO[i] = 0;
        }
        for (int i = -iSTO; i < N; i++) {
            y_STO[i] = y[i + iSTO];
        }
    }
    
    return y_STO;
}


int main (){

    int nSTO = 9;

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
    vector<complex<double>> y_CFO(y_real.size());
    for(int i=0; i<y_real.size(); i++) {
        y_CFO[i] = complex<double>(y_real[i], 0);
    }

    // Close input file
    inputFile.close();

    
    cout << "import the number of y_real :" << y_real.size() << endl;
    
    std::vector<std::complex<double>> y_CFO_STO;

    y_CFO_STO= add_STO(y_CFO,nSTO); 

    for (int n = 0; n < y_CFO_STO.size(); ++n) {
        std::cout << y_CFO_STO[n] << " ";
    }
    std::cout << std::endl;
    std::cout << "output, the number of y_CFO_STO , Output size y_CFO_STO : " << y_CFO_STO.size() << std::endl;

    return 0;
}





