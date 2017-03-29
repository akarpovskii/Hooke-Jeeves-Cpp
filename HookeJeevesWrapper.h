//
// Created by moonlightnvkz on 20.03.17.
//

#pragma once


#include <vector>
#include <functional>

namespace HookeJeevesWrapper {
    using Func = std::function<long double(const std::vector<long double>&)>;

    long double best_nearby(std::vector<long double> &delta,
                       std::vector<long double> &point, const long double prevbest,
                       const Func &f);

    size_t hooke(std::vector<long double> startpt, std::vector<long double> &endpt,
              const long double rho, const long double epsilon, const size_t itermax, const Func &f);
};


