#include <iostream>
#include "HookeJeevesWrapper.h"

double Rosenbrocks(std::vector<double> x) {
    double a, b, c;
    a = x.at(0);
    b = x.at(1);
    c = 100.0 * (b - (a * a)) * (b - (a * a));
    return (c + (1.0 - a) * (1.0 - a));
}

int main() {
    std::vector<double> startpt{-1.2, 1.0};
    std::vector<double> endpt(2);
    size_t itermax = 5000;
    double rho = 0.5;
    double epsilon = 1E-6;
    HookeJeevesWrapper::Func f (Rosenbrocks);
    std::cout << HookeJeevesWrapper::hooke(startpt, endpt, rho, epsilon, itermax, f) << std::endl;

    for (auto i : endpt) {
        std::cout << i << " ";
    }
    return 0;
}