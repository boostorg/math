// Copyright 2026 Matt Borland
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt
//
// Smoke test for the boost.math module: exercises at least one entity from
// every component exported by module/math.cxx.

#include <boost/core/lightweight_test.hpp>

// Import the STL if we can
#if defined(__cpp_lib_modules) && __cpp_lib_modules >= 202207L
import std;
#else
#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <limits>
#endif

import boost.math;

int main()
{
    // constants
    BOOST_TEST_GT(boost::math::constants::pi<double>(), 3.14);
    BOOST_TEST_LT(boost::math::constants::pi<double>(), 3.15);
    BOOST_TEST_EQ(boost::math::double_constants::half, 0.5);
    static_assert(boost::math::float_constants::half == 0.5F, "float_constants");

    // policies
    const auto pol {boost::math::policies::make_policy(
        boost::math::policies::domain_error<boost::math::policies::errno_on_error>())};
    static_cast<void>(pol);
    static_assert(boost::math::policies::is_policy<decltype(pol)>::value, "is_policy");
    BOOST_TEST_GT((boost::math::policies::digits<double, boost::math::policies::policy<>>()), 50);

    // error handling types are visible through the import
    static_assert(boost::math::policies::is_noexcept_error_policy<
        boost::math::policies::policy<>>::value == false, "default policy throws");

    // special functions
    BOOST_TEST_EQ(boost::math::tgamma(4.0), 6.0);
    BOOST_TEST_EQ(boost::math::sign(-3.5), -1);
    BOOST_TEST(boost::math::isnan(std::sqrt(-1.0)));
    BOOST_TEST_EQ(boost::math::round(2.5), 3.0);
    BOOST_TEST_EQ(boost::math::factorial<double>(5), 120.0);
    BOOST_TEST_EQ(boost::math::fpclassify(1.0), FP_NORMAL);
    {
        const auto erf_value {boost::math::erf(1.0)};
        BOOST_TEST_GT(erf_value, 0.84);
        BOOST_TEST_LT(erf_value, 0.85);
        const auto bessel_value {boost::math::cyl_bessel_j(0, 1.0)};
        BOOST_TEST_GT(bessel_value, 0.76);
        BOOST_TEST_LT(bessel_value, 0.77);
    }

    // ccmath
    static_assert(boost::math::ccmath::abs(-1) == 1, "ccmath abs");
    static_assert(boost::math::ccmath::isnan(0.0) == false, "ccmath isnan");
    static_assert(boost::math::ccmath::sqrt(4.0) == 2.0, "ccmath sqrt");

    // tools
    {
        const double coeffs[] {1.0, 2.0, 3.0};
        BOOST_TEST_EQ(boost::math::tools::evaluate_polynomial(coeffs, 2.0), 17.0);
        const auto roots {boost::math::tools::cubic_roots(1.0, -6.0, 11.0, -6.0)};
        BOOST_TEST_EQ(roots[0], 1.0);
        BOOST_TEST_EQ(roots[2], 3.0);
        const boost::math::tools::polynomial<double> p {{1.0, 1.0}};
        const auto p2 {p * p};
        BOOST_TEST_EQ(p2[1], 2.0);
    }

    // quadrature and interpolators
    {
        const auto integral {boost::math::quadrature::trapezoidal(
            [](double x) { return x * x; }, 0.0, 1.0)};
        BOOST_TEST_GT(integral, 0.333);
        BOOST_TEST_LT(integral, 0.334);
        std::vector<double> v {0.0, 1.0, 4.0, 9.0, 16.0};
        boost::math::interpolators::cardinal_cubic_b_spline<double> spline {v.data(), v.size(), 0.0, 1.0};
        const auto interp {spline(2.0)};
        BOOST_TEST_GT(interp, 3.9);
        BOOST_TEST_LT(interp, 4.1);
    }

    // statistics and optimization
    {
        std::vector<double> data {1.0, 2.0, 3.0, 4.0, 5.0};
        BOOST_TEST_EQ(boost::math::statistics::mean(data), 3.0);
        BOOST_TEST_EQ(boost::math::statistics::median(data), 3.0);
        boost::math::optimization::random_search_parameters<std::vector<double>> params {};
        params.lower_bounds = {-1.0};
        params.upper_bounds = {1.0};
        static_cast<void>(params);
    }

    // distributions
    {
        const boost::math::normal_distribution<> dist {};
        BOOST_TEST_EQ(boost::math::cdf(dist, 0.0), 0.5);
        BOOST_TEST_EQ(boost::math::median(dist), 0.0);
        BOOST_TEST_EQ(boost::math::cdf(boost::math::complement(dist, 0.0)), 0.5);
        const boost::math::students_t st {5.0};
        BOOST_TEST_GT(boost::math::pdf(st, 0.0), 0.37);
        const auto q {boost::math::quantile(boost::math::binomial(50, 0.5), 0.5)};
        BOOST_TEST_GT(q, 24.0);
        BOOST_TEST_LT(q, 26.0);
    }

    // differentiation
    {
        const auto derivative {boost::math::differentiation::finite_difference_derivative(
            [](double x) { return x * x; }, 1.0)};
        BOOST_TEST_GT(derivative, 1.99);
        BOOST_TEST_LT(derivative, 2.01);
        const auto x {boost::math::differentiation::make_fvar<double, 2>(3.0)};
        const auto y {x * x};
        BOOST_TEST_EQ(y.derivative(1), 6.0);
        // std::numeric_limits partial specialization is reachable through the import
        using fvar_type = decltype(x);
        static_assert(std::numeric_limits<fvar_type>::is_specialized, "autodiff limits");
    }

    // quaternion, octonion and complex inverse trig
    {
        const boost::math::quaternion<double> q {1.0, 2.0, 3.0, 4.0};
        BOOST_TEST_EQ(boost::math::norm(q), 30.0);
        const boost::math::octonion<double> o {1.0, 1.0};
        BOOST_TEST_GT(boost::math::abs(o), 1.41);
        const std::complex<double> z {0.5, 0.0};
        const auto az {boost::math::asin(z)};
        BOOST_TEST_GT(az.real(), 0.52);
        BOOST_TEST_LT(az.real(), 0.53);
    }

    return boost::report_errors();
}
