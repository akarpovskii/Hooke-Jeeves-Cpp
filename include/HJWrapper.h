#pragma once
#include <vector>
#include <functional>

class HJWrapper final {
public:
    using ldouble = long double;

    using Func = std::function<ldouble(const std::vector<ldouble>&)>;

    /// @param rho
    /// This is a user-supplied convergence parameter, which should be set to a
    /// value between 0.0 and 1.0. Larger values of rho give greater probability of
    /// convergence on highly nonlinear functions, at a cost of more function
    /// evaluations. Smaller values of rho reduces the number of evaluations
    /// (and the program running time), but increases the risk of non convergence.
    ///
    /// @param epsilon
    /// This is the criterion for halting the search for a minimum. When the
    /// algorithm begins to make less and less progress on each iteration, it checks
    /// the halting criterion: if the stepsize is below epsilon, terminate the
    /// iteration and return the current best estimate of the minimum.
    /// Larger values of epsilon (such as 1.0e-4) give quicker running time, but a
    /// less accurate estimate of the minimum. Smaller values of epsilon (such as
    /// 1.0e-7) give longer running time, but a more accurate estimate of the minimum.
    ///
    /// @param itermax
    /// Maximum number of iterations. A second, rarely used, halting criterion.
    /// If the algorithm uses >= itermax iterations, halt.
    HJWrapper(ldouble rho, ldouble epsilon, unsigned itermax, const Func &f);

    /// Doesn't store the end point to member field implicitly
    /// @return number of iterations performed
    unsigned hooke(std::vector<ldouble> startpt, std::vector<ldouble> &endpt) const;

    /// Stores the end point in the member field
    /// @return number of iterations performed
    unsigned hooke(std::vector<ldouble> startpt);

    /// Uses previous result as the start point and stores the result in the end point
    /// The previous result is used only if it was received using this or the function above
    /// @return number of iterations performed
    unsigned hooke();

    ldouble rho;

    ldouble epsilon;

    unsigned itermax;

    Func f;

    std::vector<ldouble> endpt;

protected:
#ifndef NDEBUG
    bool result_has_been_stored;	// Only for assert usage
#endif
    ldouble best_nearby(std::vector<ldouble> &delta,
                            std::vector<ldouble> &point, ldouble prevbest,
                            const Func &f) const;

};