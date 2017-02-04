//
// Created by moonlightnvkz on 28.01.17.
//

#pragma once

#include <vector>

typedef double (*Func)(std::vector<double> x);


typedef class HookeJeevesWrapper {
public:
    HookeJeevesWrapper() = delete;
    ~HookeJeevesWrapper() = delete;

    static unsigned find_min(Func func, std::vector<double> startpt,
                                                 std::vector<double> endpt,
                                                 double rho, double epsilon, unsigned itermax);

#if Debug
    static bool test_ax_b2_c();
    static bool test_ax_b3_c();
    static bool test_exp();
    static bool test_log();
    static bool test_sqrt();
    static bool test_polinomial();
    static bool test_rosenbrocks();
    static bool test_woods();
#endif

private:
    static double best_nearby(Func func, std::vector<double> delta,
                              std::vector<double> point, double prevbest);
} HJWrapper;

