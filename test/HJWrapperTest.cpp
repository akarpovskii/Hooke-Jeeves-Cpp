#include <gtest/gtest.h>
#include <cmath>
#include "../include/HJWrapper.h"

using ldouble = HJWrapper::ldouble;

using Func = HJWrapper::Func;

static unsigned itermax = 10000;

static ldouble rho = 0.5L;

static ldouble epsilon = 1E-6L;

TEST(HJWrapperTest, Rosenbrocks) {
    // Minimum at (1, 1)
    Func func = [] (const std::vector<ldouble> &x) {
        ldouble a, b, c;
        a = x.at(0);
        b = x.at(1);
        c = 100.0 * (b - (a * a)) * (b - (a * a));
        return (c + (1.0 - a) * (1.0 - a));
    };
    HJWrapper hj(rho, epsilon, itermax, func);
    hj.hooke({1, 1});
    ASSERT_EQ(hj.endpt.size(), 2);
    ldouble value = func(hj.endpt);
    EXPECT_NEAR(hj.endpt[0], 1.0, epsilon) << "Function value: " << value;
    EXPECT_NEAR(hj.endpt[1], 1.0, epsilon) << "Function value: " << value;
}

TEST(HJWrapperTest, SphereFunction) {
    // Minimum at (0, 0)
    Func func = [] (const std::vector<ldouble> &x) {
        return x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + x[4] * x[4];
    };
    HJWrapper hj(rho, epsilon, itermax, func);
    auto iters = hj.hooke({4, 5, 1, -2, 0});
    ASSERT_EQ(hj.endpt.size(), 5);
    ldouble value = func(hj.endpt);
    EXPECT_NEAR(hj.endpt[0], 0.0, epsilon) << "Function value: " << value;
    EXPECT_NEAR(hj.endpt[1], 0.0, epsilon) << "Function value: " << value;
    EXPECT_NEAR(hj.endpt[2], 0.0, epsilon) << "Function value: " << value;
    EXPECT_NEAR(hj.endpt[3], 0.0, epsilon) << "Function value: " << value;
    EXPECT_NEAR(hj.endpt[4], 0.0, epsilon) << "Function value: " << value;
}

TEST(HJWrapperTest, Beales) {
    // Minimum at (0, 0)
    Func func = [] (const std::vector<ldouble> &x) {
        ldouble a = x[0];
        ldouble b = x[1];
        return (1.5 - a + a*b)*(1.5 - a + a*b) +
               (2.25 - a + a*b*b)*(2.25 - a + a*b*b) +
               (2.625 - a + a*b*b*b)*(2.625 - a + a*b*b*b);
    };
    HJWrapper hj(rho, epsilon, itermax, func);
    auto iters = hj.hooke({1, 1});
    ASSERT_EQ(hj.endpt.size(), 2);
    ldouble value = func(hj.endpt);
    EXPECT_NEAR(hj.endpt[0], 3.0, epsilon) << "Function value: " << value;
    EXPECT_NEAR(hj.endpt[1], 0.5, epsilon) << "Function value: " << value;
}

ldouble plancks_law(ldouble wavelength, ldouble temperature, ldouble amplitude) {
    constexpr ldouble c = 3e8L;

    constexpr ldouble A = logl(2L * M_PIl * c);
    return A - 4 * logl(wavelength / 1000) - 1239/wavelength/((temperature+273)/11601) + amplitude;
}

TEST(HJWrapperTest, LSMPlancksLaw) {
    std::vector<ldouble> waves(500);
    constexpr const ldouble Temperature = 1042;
    constexpr const ldouble Amplitude = 38;
    std::generate(waves.begin(), waves.end(), [step = 399.5] () mutable {
        return step += 0.5;
    });
    std::vector<ldouble> spec;
    for (auto wave : waves) {
        spec.push_back(plancks_law(wave, Temperature, Amplitude));
    }

    Func func = [&] (const std::vector<ldouble> &x) {
        ldouble res = 0L;
        size_t size = waves.size();
        for (size_t i = 0; i < size; ++i) {
            ldouble s = spec.at(i);
            ldouble p = plancks_law(waves.at(i), x.at(0), x.at(1));
            res += (s - p) * (s - p);
        }
        return res;
    };

    HJWrapper hj(0.1, epsilon, itermax, func);
    auto iters = hj.hooke({1000, 35});
    ASSERT_EQ(hj.endpt.size(), 2);
    ldouble value = func(hj.endpt);
    EXPECT_NEAR(hj.endpt[0], Temperature, epsilon) << "Function value: " << value;
    EXPECT_NEAR(hj.endpt[1], Amplitude, epsilon) << "Function value: " << value;
}

TEST(HJWrapperTest, InvalidRho) {
    EXPECT_NO_THROW(HJWrapper(0.1, epsilon, itermax, {}));
    EXPECT_THROW(HJWrapper(2.0, epsilon, itermax, {}), std::invalid_argument);
}