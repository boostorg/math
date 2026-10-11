// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <cmath>
#include <random>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/type_index.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include "math_unit_test.hpp"

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

#if __has_include(<stdfloat>)
#  include <stdfloat>
#endif

using std::sqrt;
using std::abs;
using std::numeric_limits;
using boost::multiprecision::cpp_bin_float_50;

template<class Real>
void test_interpolation_condition()
{
    std::cout << "Testing interpolation condition for barycentric interpolation on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    std::mt19937 gen(4);
    boost::random::uniform_real_distribution<Real> dis(0.1f, 1);
    std::vector<Real> x(100);
    std::vector<Real> y(100);
    x[0] = dis(gen);
    y[0] = dis(gen);
    for (size_t i = 1; i < x.size(); ++i)
    {
        x[i] = x[i-1] + dis(gen);
        y[i] = dis(gen);
    }

    boost::math::interpolators::barycentric_rational<Real> interpolator(x.data(), y.data(), y.size());

    for (size_t i = 0; i < x.size(); ++i)
    {
        Real z = interpolator(x[i]);
        CHECK_ULP_CLOSE(y[i], z, 2);
    }

    // Make sure that the move constructor does the same thing:
    std::vector<Real> x_copy = x;
    std::vector<Real> y_copy = y;
    boost::math::interpolators::barycentric_rational<Real> move_interpolator(std::move(x), std::move(y));

    for (size_t i = 0; i < x_copy.size(); ++i)
    {
        Real z = move_interpolator(x_copy[i]);
        CHECK_ULP_CLOSE(y_copy[i], z, 2);
    }
}

template<class Real>
void test_interpolation_condition_high_order()
{
    std::cout << "Testing interpolation condition in high order for barycentric interpolation on type " << boost::typeindex::type_id<Real>().pretty_name()  << "\n";
    std::mt19937 gen(5);
    boost::random::uniform_real_distribution<Real> dis(0.1f, 1);
    std::vector<Real> x(100);
    std::vector<Real> y(100);
    x[0] = dis(gen);
    y[0] = dis(gen);
    for (size_t i = 1; i < x.size(); ++i)
    {
        x[i] = x[i-1] + dis(gen);
        y[i] = dis(gen);
    }

    // Order 5 approximation:
    boost::math::interpolators::barycentric_rational<Real> interpolator(x.data(), y.data(), y.size(), 5);

    for (size_t i = 0; i < x.size(); ++i)
    {
        Real z = interpolator(x[i]);
        CHECK_ULP_CLOSE(y[i], z, 2);
    }
}


template<class Real>
void test_constant()
{
    std::cout << "Testing that constants are interpolated correctly using barycentric interpolation on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    std::mt19937 gen(6);
    boost::random::uniform_real_distribution<Real> dis(0.1f, 1);
    std::vector<Real> x(100);
    std::vector<Real> y(100);
    Real constant = -8;
    x[0] = dis(gen);
    y[0] = constant;
    for (size_t i = 1; i < x.size(); ++i)
    {
        x[i] = x[i-1] + dis(gen);
        y[i] = y[0];
    }

    boost::math::interpolators::barycentric_rational<Real> interpolator(x.data(), y.data(), y.size());

    for (size_t i = 0; i < x.size(); ++i)
    {
        // Don't evaluate the constant at x[i]; that's already tested in the interpolation condition test.
        Real t = x[i] + dis(gen);
        Real z = interpolator(t);
        CHECK_ABSOLUTE_ERROR(constant, z, Real((sqrt(numeric_limits<Real>::epsilon()))*abs(constant)));
        CHECK_ABSOLUTE_ERROR(0, interpolator.prime(t), sqrt(numeric_limits<Real>::epsilon()));
    }
}

template<class Real>
void test_constant_high_order()
{
    std::cout << "Testing that constants are interpolated correctly in high order using barycentric interpolation on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    std::mt19937 gen(7);
    boost::random::uniform_real_distribution<Real> dis(0.1f, 1);
    std::vector<Real> x(100);
    std::vector<Real> y(100);
    Real constant = 5;
    x[0] = dis(gen);
    y[0] = constant;
    for (size_t i = 1; i < x.size(); ++i)
    {
        x[i] = x[i-1] + dis(gen);
        y[i] = y[0];
    }

    // Set interpolation order to 7:
    boost::math::interpolators::barycentric_rational<Real> interpolator(x.data(), y.data(), y.size(), 7);

    for (size_t i = 0; i < x.size(); ++i)
    {
        Real t = x[i] + dis(gen);
        Real z = interpolator(t);
        CHECK_ABSOLUTE_ERROR(constant, z, Real(10*(sqrt(numeric_limits<Real>::epsilon()))*abs(constant)));
        CHECK_ABSOLUTE_ERROR(0, interpolator.prime(t), 100*sqrt(numeric_limits<Real>::epsilon()));
    }
}


template<class Real>
void test_runge()
{
    std::cout << "Testing interpolation of Runge's 1/(1+25x^2) function using barycentric interpolation on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    std::mt19937 gen(8);
    boost::random::uniform_real_distribution<Real> dis(0.005f, 0.01f);
    std::vector<Real> x(100);
    std::vector<Real> y(100);
    x[0] = -2;
    y[0] = 1/(1+25*x[0]*x[0]);
    for (size_t i = 1; i < x.size(); ++i)
    {
        x[i] = x[i-1] + dis(gen);
        y[i] = 1/(1+25*x[i]*x[i]);
    }

    boost::math::interpolators::barycentric_rational<Real> interpolator(x.data(), y.data(), y.size(), 5);

    for (size_t i = 0; i < x.size(); ++i)
    {
        Real t = x[i];
        Real z = interpolator(t);
        CHECK_ABSOLUTE_ERROR(y[i], z, Real(0.0003*abs(y[i])));
        Real z_prime = interpolator.prime(t);
        Real num = -50*t;
        Real denom = (1+25*t*t)*(1+25*t*t);
        if (abs(num/denom) > 0.00001)
        {
            CHECK_ABSOLUTE_ERROR(num/denom, z_prime, Real(0.03*abs(num/denom)));
        }
    }


    Real tol = 0.0001;
    for (size_t i = 0; i < x.size(); ++i)
    {
        Real t = x[i] + dis(gen);
        Real z = interpolator(t);
        CHECK_ABSOLUTE_ERROR(1/(1+25*t*t), z, Real(tol/100*abs(1/(1+25*t*t))));
        Real z_prime = interpolator.prime(t);
        Real num = -50*t;
        Real denom = (1+25*t*t)*(1+25*t*t);
        Real runge_prime = num/denom;

        if (abs(runge_prime) > 0 && abs(z_prime - runge_prime)/abs(runge_prime) > tol)
        {
            std::cout << "Error too high for t = " << t << " which is a distance " << t - x[i] << " from node " << i << "/" << x.size() << " associated with data (" << x[i] << ", " << y[i] << ")\n";
            CHECK_ABSOLUTE_ERROR(runge_prime, z_prime, Real(tol*abs(runge_prime)));
        }
    }
}

template<class Real>
void test_weights()
{
    std::cout << "Testing weights are calculated correctly using barycentric interpolation on type " << boost::typeindex::type_id<Real>().pretty_name() << "\n";

    std::mt19937 gen(9);
    boost::random::uniform_real_distribution<Real> dis(0.005, 0.01);
    std::vector<Real> x(100);
    std::vector<Real> y(100);
    x[0] = -2;
    y[0] = 1/(1+25*x[0]*x[0]);
    for (size_t i = 1; i < x.size(); ++i)
    {
        x[i] = x[i-1] + dis(gen);
        y[i] = 1/(1+25*x[i]*x[i]);
    }

    boost::math::interpolators::detail::barycentric_rational_imp<Real> interpolator(x.data(), x.data() + x.size(), y.data(), 0);

    for (size_t i = 0; i < x.size(); ++i)
    {
        Real w = interpolator.weight(i);
        if (i % 2 == 0)
        {
            CHECK_ABSOLUTE_ERROR(1, w, Real(1e-07));
        }
        else
        {
            CHECK_ABSOLUTE_ERROR(-1, w, Real(1e-07));
        }
    }

    // d = 1:
    interpolator = boost::math::interpolators::detail::barycentric_rational_imp<Real>(x.data(), x.data() + x.size(), y.data(), 1);

    for (size_t i = 1; i < x.size() -1; ++i)
    {
        Real w = interpolator.weight(i);
        Real w_expect = 1/(x[i] - x[i - 1]) + 1/(x[i+1] - x[i]);
        if (i % 2 == 0)
        {
            CHECK_ABSOLUTE_ERROR(-w_expect, w, Real(1e-07*abs(-w_expect)));
        }
        else
        {
            CHECK_ABSOLUTE_ERROR(w_expect, w, Real(1e-07*abs(w_expect)));
        }
    }

}


void barycentric_rational()
{
    // The tests took too long at the higher precisions.
    // They still pass, but the CI system is starting to time out,
    // so I figured it'd be polite to comment out the most expensive tests.
    
    #ifdef __STDCPP_FLOAT32_T__
    
    test_constant<std::float32_t>();
    //test_constant_high_order<std::float32_t>();
    test_interpolation_condition<std::float32_t>();
    //test_interpolation_condition_high_order<std::float32_t>();
    
    #else
    
    test_constant<float>();
    //test_constant_high_order<float>();
    test_interpolation_condition<float>();
    //test_interpolation_condition_high_order<float>();
    
    #endif
    
    #ifdef __STDCPP_FLOAT64_T__
    
    test_weights<std::float64_t>();
    //test_constant<std::float64_t>();
    test_constant_high_order<std::float64_t>();
    test_interpolation_condition<std::float64_t>();
    test_interpolation_condition_high_order<std::float64_t>();
    test_runge<std::float64_t>();
    
    #else
    
    test_weights<double>();
    //test_constant<double>();
    test_constant_high_order<double>();
    test_interpolation_condition<double>();
    test_interpolation_condition_high_order<double>();
    test_runge<double>();
    
    #endif

    test_constant<long double>();
    //test_constant_high_order<long double>();
    //test_interpolation_condition<long double>();
    //test_interpolation_condition_high_order<long double>();
    //test_runge<long double>();

    //test_constant<cpp_bin_float_50>();
    //test_constant_high_order<cpp_bin_float_50>();
    //test_interpolation_condition<cpp_bin_float_50>();
    //test_interpolation_condition_high_order<cpp_bin_float_50>();
    //test_runge<cpp_bin_float_50>();

#ifdef BOOST_HAS_FLOAT128
    //test_interpolation_condition<boost::multiprecision::float128>();
    //test_constant<boost::multiprecision::float128>();
    //test_constant_high_order<boost::multiprecision::float128>();
    //test_interpolation_condition_high_order<boost::multiprecision::float128>();
    //test_runge<boost::multiprecision::float128>();
#endif

}

int main()
{
    barycentric_rational();
    return boost::math::test::report_errors();
}
