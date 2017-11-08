#include <cassert>
#include <cmath>
#include <utility>
#include "../include/HJWrapper.h"


HJWrapper::HJWrapper(ldouble rho, ldouble epsilon, unsigned itermax, const Func & f) {
    if (fabsl(rho) > 1.0l) {
        throw std::invalid_argument("rho must be -1 < rho < 1");
    }
    this->rho = rho;
    this->epsilon = epsilon;
    this->itermax = itermax;
    this->f = f;
}

HJWrapper::ldouble HJWrapper::best_nearby(std::vector<ldouble> &delta, std::vector<ldouble> &point,
                                            ldouble prevbest, const Func &f) const {
    size_t size = point.size();
    assert(delta.size() == point.size());

    std::vector<ldouble> tmp(point.begin(), point.end());
    ldouble fmin = prevbest;

    for (size_t i = 0; i < size; ++i) {
        tmp.at(i) = point.at(i) + delta.at(i);
        ldouble newf = f(tmp);
        if (newf < fmin) {
            fmin = newf;
        } else {
            delta.at(i) *= -1;
            tmp.at(i) = point.at(i) + delta.at(i);
            newf = f(tmp);
            if (newf < fmin) {
                fmin = newf;
            } else {
                tmp.at(i) = point.at(i);
            }
        }
    }
    point = std::move(tmp);
    return fmin;
}

unsigned HJWrapper::hooke(const std::vector<ldouble> &startpt_, std::vector<ldouble> &endpt) const {
    auto startpt = startpt_;
    size_t size = startpt.size();
    std::vector<ldouble> delta(size);
    unsigned iterations = 0;

    for (size_t i = 0; i < size; ++i) {
        delta.at(i) = fabsl(startpt.at(i) * rho);
        if (delta.at(i) < epsilon) {
            delta.at(i) = rho;
        }
    }

    ldouble step_length = rho;
    ldouble fbefore = f(startpt);
    while (iterations < itermax && step_length > epsilon) {
        ++iterations;
        /* find best new point, one coord at a time */
        endpt = startpt;
        ldouble newf = best_nearby(delta, endpt, fbefore, f);

        /* if we made some improvements, pursue that direction */
        bool keep = true;
        while (newf < fbefore && keep) {
            for (size_t i = 0; i < size; ++i) {
                /* firstly, arrange the sign of delta[] */
                if (endpt.at(i) < startpt.at(i)) {
                    delta.at(i) = -fabsl(delta.at(i));
                } else {
                    delta.at(i) = fabsl(delta.at(i));
                }
                /* now, move further in this direction */
                ldouble dx = endpt.at(i) - startpt.at(i);
                startpt.at(i) = endpt.at(i);
                endpt.at(i) = endpt.at(i) + dx;
            }
            fbefore = newf;
            newf = best_nearby(delta, endpt, fbefore, f);
            /* if the further (optimistic) move was bad.... */
            if (newf > fbefore)
                break;
            /* make sure that the differences between the new */
            /* and the old points are due to actual */
            /* displacements; beware of roundoff errors that */
            /* might cause newf < fbefore */
            for (unsigned i = 0; i < size; ++i) {
                keep = true;
                if (fabsl(endpt.at(i) - startpt.at(i)) > 0.5L * fabsl(delta.at(i))) {
                    break;
                } else {
                    keep = false;
                }
            }
        }
        if (step_length > epsilon && newf > fbefore) {
            step_length *= rho;
            for (auto &d : delta) {
                d *= rho;
            }
        }
    }
    return iterations;
}

unsigned HJWrapper::hooke(const std::vector<ldouble> &startpt) {
    unsigned iters = hooke(startpt, endpt);
#ifndef NDEBUG
    result_has_been_stored = true;
#endif
    return iters;
}

unsigned HJWrapper::hooke() {
    assert(result_has_been_stored);
    return hooke(endpt, endpt);
}
