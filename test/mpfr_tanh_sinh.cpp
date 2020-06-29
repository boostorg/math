/*
 * Copyright Nick Thompson, 2020
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/mpfr.hpp>

using boost::math::quadrature::tanh_sinh;
using boost::math::constants::pi;
using boost::multiprecision::mpfr_float;

template <class Real>
const tanh_sinh<Real>& get_integrator()
{
   static const tanh_sinh<Real> integrator(15);
   return integrator;
}

template<class Real>
void test_linear()
{
    using std::sqrt;
    std::cout << "Testing linear functions are integrated properly by tanh_sinh on type " << boost::core::demangle(typeid(Real).name()) << "\n";
    auto integrator = get_integrator<Real>();
    auto f = [](const Real& x)->Real
    {
       return 5*x + 7;
    };
    Real error;
    Real L1;
    Real Q = integrator.integrate(f, (Real) 0, (Real) 1, sqrt(std::numeric_limits<Real>::epsilon()), &error, &L1);
    Real expected = 9.5;
    CHECK_ULP_CLOSE(expected, Q, 10);
    CHECK_ULP_CLOSE(expected, L1, 10);
    Q = integrator.integrate(f, (Real) 1, (Real) 0,sqrt(std::numeric_limits<Real>::epsilon()), &error, &L1);
    Real negated = -expected;
    CHECK_ULP_CLOSE(negated, Q, 10);
    CHECK_ULP_CLOSE(expected, L1, 10);
    Q = integrator.integrate(f, (Real) 1, (Real) 1, sqrt(std::numeric_limits<Real>::epsilon()), &error, &L1);
    CHECK_EQUAL(Real(0), Q);
}


int main()
{
    test_linear<mpfr_float>();
    return boost::math::test::report_errors();
}
