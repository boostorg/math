// Copyright John Maddock 2006.
// Copyright Paul A. Bristow 2007, 2009
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error

#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/tools/stats.hpp>
#include <boost/math/tools/test.hpp>
#include <boost/math/tools/big_constant.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include "functor.hpp"

#include "handle_test_result.hpp"
#include "table_type.hpp"

#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>

template <class Real, class T>
void do_test_2F0(const T& data, const char* type_name, const char* test_name)
{
   typedef Real                   value_type;

   typedef value_type(*pg)(value_type, value_type, value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
   pg funcp = boost::math::hypergeometric_0F1<value_type, value_type>;
#else
   pg funcp = boost::math::hypergeometric_1F1;
#endif

   boost::math::tools::test_result<value_type> result;

   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   //
   // test hypergeometric_2F0 against data:
   //
   result = boost::math::tools::test_hetero<Real>(
      data,
      bind_func<Real>(funcp, 0, 1, 2),
      extract_result<Real>(3));
   handle_test_result(result, data[result.worst()], result.worst(), type_name, "hypergeometric_1F1", test_name);
   std::cout << std::endl;
}

#ifndef SC_
#define SC_(x) BOOST_MATH_BIG_CONSTANT(T, 1000000, x)
#endif

template <class T>
void test_spots1(T, const char* type_name)
{
#include "hypergeometric_1F1.ipp"

   do_test_2F0<T>(hypergeometric_1F1, type_name, "Integer a values");

#include "hypergeometric_1F1_small_random.ipp"

   do_test_2F0<T>(hypergeometric_1F1_small_random, type_name, "Small random values");
}

template <class T>
void test_spots2(T, const char* type_name)
{
#include "hypergeometric_1F1_big.ipp"

   do_test_2F0<T>(hypergeometric_1F1_big, type_name, "Large random values");
}

template <class T>
void test_spots(T z, const char* type_name)
{
   test_spots1(z, type_name);
   test_spots2(z, type_name);
}


// Tests the Mellin transform formula given here: https://dlmf.nist.gov/13.10, Equation 13.10.10
template <class Real>
void test_hypergeometric_mellin_transform()
{
    using boost::math::hypergeometric_1F1;
    using boost::math::quadrature::exp_sinh;
    using boost::math::tgamma;
    using std::pow;

    // Constraint: 0 < lambda < a.
    Real lambda = 0.5;
    Real a = 1;
    Real b = 3;
    auto f = [&](Real t)->Real { return pow(t, lambda - 1)*hypergeometric_1F1(a, b, -t); };

    auto integrator = exp_sinh<double>();
    Real computed = integrator.integrate(f, boost::math::tools::epsilon<Real>());
    Real expected = tgamma(b)*tgamma(lambda)*tgamma(a-lambda)/(tgamma(a)*tgamma(b-lambda));

    Real tol = boost::math::tools::epsilon<Real>() * 10;
    BOOST_CHECK_CLOSE(computed, expected, tol);
}


// Tests the Laplace transform formula given here: https://dlmf.nist.gov/13.10, Equation 13.10.4
template <class Real>
void test_hypergeometric_laplace_transform()
{
    using boost::math::hypergeometric_1F1;
    using boost::math::quadrature::exp_sinh;
    using boost::math::tgamma;
    using std::pow;
    using std::exp;

    // Set a = 1 blows up for some reason . . .
    Real a = -1;
    Real b = 3;
    Real z = 1.5;
    auto f = [&](Real t)->Real { return exp(-z*t)*pow(t, b - 1)*hypergeometric_1F1(a, b, t); };

    auto integrator = exp_sinh<double>();
    Real computed = integrator.integrate(f, boost::math::tools::epsilon<Real>());
    Real expected = tgamma(b)/(pow(z,b)*pow(1-1/z, a));

    Real tol = boost::math::tools::epsilon<Real>() * 200;
    BOOST_CHECK_CLOSE(computed, expected, tol);
}
