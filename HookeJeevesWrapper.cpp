//
// Created by moonlightnvkz on 20.03.17.
//

#include <cmath>
#include "HookeJeevesWrapper.h"

long double HookeJeevesWrapper::best_nearby(std::vector<long double> &delta, std::vector<long double> &point,
                                       const long double prevbest, const Func &f) {
    if (delta.size() != point.size()) {
        throw std::out_of_range("Sizes of delta and point is different");
    }
    std::vector<long double> tmp(point.begin(), point.end());
    long double fmin = prevbest, newf;

    const size_t size = point.size();
    for (size_t i = 0; i < size; i++) {
        tmp.at(i) = point.at(i) + delta.at(i);
        newf = f(tmp);
        if (newf < fmin)
            fmin = newf;
        else {
            delta.at(i) *= -1;
            tmp.at(i) = point.at(i) + delta.at(i);
            newf = f(tmp);
            if (newf < fmin)
                fmin = newf;
            else
                tmp.at(i) = point.at(i);
        }
    }
    point = std::move(tmp);
    return fmin;
}

size_t HookeJeevesWrapper::hooke(std::vector<long double> startpt, std::vector<long double> &endpt,
                                 const long double rho, const long double epsilon, const size_t itermax, const Func &f) {
    size_t size = startpt.size();
    endpt.resize(size);
    std::vector<long double> delta(size);
    size_t iterations = 0;

    constexpr long double Error = 1E-12;
    for (size_t i = 0; i < size; i++) {
        delta.at(i) = std::fabs(startpt.at(i) * rho);
        if (delta.at(i) < Error)
            delta.at(i) = rho;
    }

    long double step_length = rho;
    long double fbefore = f(startpt), newf;
    bool keep;
    while (iterations < itermax && step_length > epsilon) {
        iterations++;
        /* find best new point, one coord at a time */
        std::copy(startpt.begin(), startpt.end(), endpt.begin());
        newf = best_nearby(delta, endpt, fbefore, f);

        /* if we made some improvements, pursue that direction */
        keep = true;
        while (newf < fbefore && keep) {
            for (size_t i = 0; i < size; i++) {
                /* firstly, arrange the sign of delta[] */
                if (endpt.at(i) <= startpt.at(i))
                    delta.at(i) = -fabsl(delta.at(i));
                else
                    delta.at(i) = fabsl(delta.at(i));
                /* now, move further in this direction */
                long double dx = endpt.at(i) - startpt.at(i);
                startpt.at(i) = endpt.at(i);
                endpt.at(i) = endpt.at(i) + dx;
            }
            fbefore = newf;
            newf = best_nearby(delta, endpt, fbefore, f);
            /* if the further (optimistic) move was bad.... */
            if (newf >= fbefore)
                break;
            /* make sure that the differences between the new */
            /* and the old points are due to actual */
            /* displacements; beware of roundoff errors that */
            /* might cause newf < fbefore */
            for (size_t i = 0; i < size; i++) {
                keep = true;
                if (fabsl(endpt.at(i) - startpt.at(i)) > 0.5L * fabsl(delta.at(i)))
                    break;
                else
                    keep = false;
            }
        }
        if (step_length >= epsilon && newf >= fbefore) {
            step_length *= rho;
            for (auto &d : delta) {
                d *= rho;
            }
        }
    }
    endpt = std::move(startpt);
    return iterations;
}
