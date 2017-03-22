//
// Created by moonlightnvkz on 20.03.17.
//

#pragma once


#include <vector>
#include <functional>

namespace HookeJeevesWrapper {
    using Func = std::function<double(std::vector<double>)>;

    double best_nearby(std::vector<double> &delta,
                       std::vector<double> &point, const double prevbest,
                       const Func &f);

    size_t hooke(std::vector<double> startpt, std::vector<double> &endpt,
              double rho, double epsilon, size_t itermax, const Func &f);
};


