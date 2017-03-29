#include <iostream>
#include <cmath>
#include <fstream>
#include "HookeJeevesWrapper.h"

long double Rosenbrocks(std::vector<long double> x) {
    long double a, b, c;
    a = x.at(0);
    b = x.at(1);
    c = 100.0 * (b - (a * a)) * (b - (a * a));
    return (c + (1.0 - a) * (1.0 - a));
}

long double sinus(std::vector<long double> x) {
    return sinl(x[0]);
}

long double x2(std::vector<long double> x) {
    return (x[0] - 2) * (x[0] - 2);
}

long double x3(std::vector<long double> x) {
    return (x[0] - 1) * (x[0] - 2) * (x[0] - 3);
}

long double x2_2d(std::vector<long double> x) {
    return (x[0] - 2) * (x[0] - 2) + (x[1] - 3) * (x[1] - 3);
}

long double plancks_law(long double wavelength, long double temperature, long double amplitude) {
    constexpr long double h = 6.625e-34L;
    constexpr long double k = 1.38e-23L;
    constexpr long double c = 3e8L;

    constexpr long double A = logl(2L * M_PIl * c);
    return A - 4 * logl(wavelength / 1000) - 1239/wavelength/((temperature+273)/11601) + amplitude;
}

std::vector<long double> waves;
std::vector<long double> spec;

long double lsm_plancks_law(const std::vector<long double> &x) {
    long double res = 0L;
    size_t size = waves.size();
    for (size_t i = 240; i < size; i += 1) {
        long double s = spec.at(i);
        long double p = plancks_law(waves.at(i), x.at(0), x.at(1));
        res += (s - p) * (s - p);
    }
    return res;
}

bool read_data(std::vector<long double> &x, std::vector<long double> &y) {
    std::ifstream is("1042.dat");    // 657.dat and 1042.dat
    if (!is) {
        return false;
    }
    while (is) {
        long double _x, _y;
        is >> _x >> _y;
        x.push_back(_x);
        y.push_back(_y);
    }
    return true;
}

int main() {
    if (!read_data(waves, spec)) {
        return 1;
    }
    std::vector<long double> startpt{0, 0};
    std::vector<long double> endpt;
    size_t itermax = 10000;
    long double rho = 0.5L;
    long double epsilon = 1E-15L;
    HookeJeevesWrapper::Func f(lsm_plancks_law);
    std::cout << HookeJeevesWrapper::hooke(startpt, endpt, rho, epsilon, itermax, f) << std::endl;

    for (auto e : endpt) {
        std::cout << e << " ";
    }
    std::cout << std::endl << lsm_plancks_law(endpt);
    return 0;
}