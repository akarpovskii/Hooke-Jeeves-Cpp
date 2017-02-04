//
// Created by moonlightnvkz on 28.01.17.
//

#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <algorithm>
#include <cassert>
#include "HookeJeevesWrapper.h"

using namespace std;

unsigned HookeJeevesWrapper::find_min(Func func, vector<double> startpt,
                                      vector<double> endpt,
                                      double rho, double epsilon, unsigned itermax) {
    assert(endpt.size() == startpt.size());
    vector<double> xbefore = startpt;
    vector<double> newx = startpt;
    vector<double> delta(startpt.size());

    for (size_t i = 0; i < startpt.size(); ++i) {
        delta.at(i) = fabs(startpt.at(i) * rho);
        // Was: delta[i] == 0.0 ???
        if (delta.at(i) < epsilon) {
            delta.at(i) = rho;
        }
    }
    int iadj = 0;
    double steplength = rho;
    unsigned iters = 0;
    double fbefore = func(newx);
    double newf;

    while (iters < itermax && steplength > epsilon) {
        ++iters;
        ++iadj;
        newx = xbefore;
        newf = best_nearby(func, delta, newx, fbefore);
        /* if we made some improvements, pursue that direction */
        bool keep = true;
        while (newf < fbefore && keep) {
            iadj = 0;
            for (size_t i = 0; i < xbefore.size(); ++i) {
                /* firstly, arrange the sign of delta[] */
                if (newx.at(i) <= xbefore.at(i)) {
                    delta.at(i) = -fabs(delta.at(i));
                } else {
                    delta.at(i) = fabs(delta.at(i));
                }
                double temp = xbefore.at(i);
                xbefore.at(i) = newx.at(i);
                newx.at(i) = newx.at(i) + newx.at(i) - temp;
            }
            fbefore = newf;
            newf = best_nearby(func, delta, newx, fbefore);
            /* if the further (optimistic) move was bad.... */
            if (newf >= fbefore) {
                break;
            }
            /* make sure that the differences between the new */
            /* and the old points are due to actual */
            /* displacements; beware of roundoff errors that */
            /* might cause newf < fbefore */
            keep = false;
            for (size_t i = 0; i < xbefore.size(); ++i) {
                keep = true;
                if (fabs(newx[i] - xbefore[i]) > 0.5 * fabs(delta[i])) {
                    break;
                } else {
                    keep = false;
                }
            }
        }
        if (steplength >= epsilon && newf >= fbefore) {
            steplength *= rho;
            for (auto &&d : delta) {
                d *= rho;
            }
        }
    }
    for (size_t i = 0; i < xbefore.size(); ++i) {
        endpt[i] = xbefore[i];
    }
    return iters;
}

double HookeJeevesWrapper::best_nearby(Func func, vector<double> delta,
                                       vector<double> point, double prevbest) {
    assert(point.size() == delta.size());
    vector<double> z(delta.size());
    double minf = prevbest;
    for (size_t i = 0; i < delta.size(); ++i) {
        z.at(i) = point.at(i) + delta.at(i);
        double ftemp = func(z);
        if (ftemp < minf) {
            minf = ftemp;
        } else {
            delta.at(i) = -delta.at(i);
            z.at(i) = point.at(i) + delta.at(i);
            ftemp = func(z);
            if (ftemp < minf) {
                minf = ftemp;
            } else {
                z.at(i) = ftemp;
            }
        }
    }
    for (size_t i = 0; i < z.size(); ++i) {
        point[i] = z[i];
    }
    return minf;
}

#if Debug

static void fill_measured_values(double *x_measured, double *y_measured,
                                 int interval_length, std::function<double (double)> f) {
    for (int i = 0; i < interval_length; ++i) {
        x_measured[i] = i;
        y_measured[i] = f(x_measured[i]);
    }
}

bool HookeJeevesWrapper::test_ax_b2_c() {
    const size_t interval_length = 101;
    const double a = 13;
    const double b = 24;
    const double c = 0;

    // The only way not to capture the variables (func lambda) is to make them static
    static double x_measured[interval_length], y_measured[interval_length];
    fill_measured_values(x_measured, y_measured, interval_length, [a, b, c](double x) {
        return a * (x - b)*(x - b) + c;
    });

    auto func = [](double *coefficients, size_t n) {
        double sum = 0;
        for (size_t i = 0; i < interval_length; ++i) {
            double summand = (y_measured[i] -
                              (coefficients[0] * (x_measured[i] - coefficients[1]) * (x_measured[i] - coefficients[1]) +
                               coefficients[2]));
            sum += summand * summand;
        }
        return sum;
    };

    const size_t nvars = 3;
    double startpt[nvars] = {0};
    double endpt[nvars];
    double rho = 0.5, epsilon = 0.000001;
    unsigned itermax = 100000;

    std::cout << "Iterations required: "
              << find_min(func, startpt, endpt, nvars, rho, epsilon, itermax) << std::endl;
    std::cout << "Coefficients: " << endpt[0] << ", " << endpt[1] << ", " << endpt[2] << std::endl;
    return fabs(endpt[0] - a) < epsilon &&
           fabs(endpt[1] - b) < epsilon;
}

bool HookeJeevesWrapper::test_ax_b3_c() {
    const size_t interval_length = 101;
    const double a = 1;
    const double b = 2;
    const double c = 3;

    // The only way not to capture the variables (func lambda) is to make them static
    static double x_measured[interval_length], y_measured[interval_length];
    fill_measured_values(x_measured, y_measured, interval_length, [a, b, c](double x) {
        return a * (x - b)*(x - b)*(x - b) + c;
    });

    auto func = [](double *coefficients, size_t n) {
        double sum = 0;
        for (size_t i = 0; i < interval_length; ++i) {
            double summand = (y_measured[i] -
                              (coefficients[0] * (x_measured[i] - coefficients[1]) * (x_measured[i] - coefficients[1]) *
                               (x_measured[i] - coefficients[1]) + coefficients[2]));
            sum += summand * summand;
        }
        return sum;
    };

    const size_t nvars = 3;
    double startpt[nvars] = {1, 2, 3};
    double endpt[nvars];
    double rho = 0.7, epsilon = 0.000001;
    unsigned itermax = 100000;

    std::cout << "Iterations required: "
              << find_min(func, startpt, endpt, nvars, rho, epsilon, itermax) << std::endl;
    std::cout << "Coefficients: " << endpt[0] << ", " << endpt[1] << ", " << endpt[2] << std::endl;
    return fabs(endpt[0] - a) < epsilon &&
           fabs(endpt[1] - b) < epsilon;
}

bool HookeJeevesWrapper::test_exp() {
    const size_t interval_length = 11;
    const double a = 1;
    const double b = 1;
    const double c = 0;
    const double d = 0;

    // The only way not to capture the variables (func lambda) is to make them static
    static double x_measured[interval_length], y_measured[interval_length];
    fill_measured_values(x_measured, y_measured, interval_length, [a, b, c, d](double x) {
        return a * exp(b * x + c) + d;
    });

    auto func = [](double *coefficients, size_t n) {
        double sum = 0;
        for (size_t i = 0; i < interval_length; ++i) {
            double summand = (y_measured[i] -
                              (coefficients[0] * exp(coefficients[1] * x_measured[i] + coefficients[2]) +
                               coefficients[3]));
            sum += summand * summand;
        }
        return sum;
    };

    const size_t nvars = 4;
    double startpt[nvars] = {0};
    double endpt[nvars];
    double rho = 0.85, epsilon = 1E-6;
    unsigned itermax = 100000;

    std::cout << "Iterations required: "
              << find_min(func, startpt, endpt, nvars, rho, epsilon, itermax) << std::endl;
    std::cout << "Coefficients: " << endpt[0] << ", " << endpt[1] << ", " << endpt[2] << ", " << endpt[3] << std::endl;
    return fabs(endpt[0] - a) < epsilon &&
           fabs(endpt[1] - b) <  epsilon;
}

bool HookeJeevesWrapper::test_log() {
    return false;
}

bool HookeJeevesWrapper::test_sqrt() {
    return false;
}

bool HookeJeevesWrapper::test_polinomial() {
    return false;
}

bool HookeJeevesWrapper::test_rosenbrocks() {
    return false;
}

bool HookeJeevesWrapper::test_woods() {
    return false;
}
#endif