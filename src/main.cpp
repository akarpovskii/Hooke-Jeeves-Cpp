#include <iostream>
#include <cmath>
#include <fstream>
#include "../include/HJWrapper.h"

using HJWrapper::ldouble;

ldouble Rosenbrocks(const std::vector<ldouble> &x) {
    ldouble a, b, c;
    a = x.at(0);
    b = x.at(1);
    c = 100.0 * (b - (a * a)) * (b - (a * a));
    return (c + (1.0 - a) * (1.0 - a));
}

ldouble sinus(const std::vector<ldouble> &x) {
    return sinl(x[0]);
}

ldouble x2(const std::vector<ldouble> &x) {
    return (x[0] - 2) * (x[0] - 2);
}

ldouble x3(const std::vector<ldouble> &x) {
    return (x[0] - 1) * (x[0] - 2) * (x[0] - 3);
}

ldouble x2_2d(const std::vector<ldouble> &x) {
    return (x[0] - 2) * (x[0] - 2) + (x[1] - 3) * (x[1] - 3);
}

ldouble plancks_law(const ldouble wavelength, const ldouble temperature, const ldouble amplitude) {
    constexpr ldouble c = 3e8L;
    constexpr ldouble A = logl(2L * M_PIl * c);
    return A - 4 * logl(wavelength / 1000) - 1239/wavelength/((temperature+273)/11601) + amplitude;
}

std::vector<ldouble> waves;
std::vector<ldouble> spec;

ldouble lsm_plancks_law(const std::vector<ldouble> &x) {
    ldouble res = 0L;
    size_t size = waves.size();
    for (size_t i = 240; i < size; i += 1) {
        ldouble s = spec.at(i);
        ldouble p = plancks_law(waves.at(i), x.at(0), x.at(1));
        res += (s - p) * (s - p);
    }
    return res;
}

bool read_data(std::vector<ldouble> &x, std::vector<ldouble> &y) {
    std::ifstream is("1042.dat");    // 657.dat and 1042.dat
    if (!is) {
        return false;
    }
    while (is) {
        ldouble _x, _y;
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
    std::vector<ldouble> startpt{0, 0};
    constexpr size_t itermax = 10000;
    constexpr ldouble rho = 0.5L;
    constexpr ldouble epsilon = 1E-15L;
    HJWrapper::Func f(lsm_plancks_law);
    HJWrapper hjWrapper(rho, epsilon, itermax, f);
    std::cout << hjWrapper.hooke(startpt) << std::endl;

    for (auto e : hjWrapper.endpt) {
        std::cout << e << " ";
    }
    std::cout << std::endl << lsm_plancks_law(hjWrapper.endpt);
    return 0;
}