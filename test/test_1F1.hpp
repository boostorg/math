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

#ifdef BOOST_MSVC
#pragma warning(disable:4127)
#endif

template <class Real, class T>
void do_test_1F1(const T& data, const char* type_name, const char* test_name)
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

   do_test_1F1<T>(hypergeometric_1F1, type_name, "Integer a values");

#include "hypergeometric_1F1_small_random.ipp"

   do_test_1F1<T>(hypergeometric_1F1_small_random, type_name, "Small random values");
}

template <class T>
void test_spots2(T, const char* type_name)
{
#include "hypergeometric_1F1_big.ipp"

   do_test_1F1<T>(hypergeometric_1F1_big, type_name, "Large random values");
}

template <class T>
void test_spots3(T, const char* type_name)
{
#include "hypergeometric_1F1_big_double_limited.ipp"

   do_test_1F1<T>(hypergeometric_1F1_big_double_limited, type_name, "Large random values - double limited precision");
}

template <class T>
void test_spots4(T, const char* type_name)
{
#include "hypergeometric_1F1_big_unsolved.ipp"

   do_test_1F1<T>(hypergeometric_1F1_big, type_name, "Large random values - unsolved domains");
}

template <class T>
void test_spots5(T, const char* type_name)
{
   std::cout << "Testing special cases for type " << type_name << std::endl;
   BOOST_MATH_STD_USING
   //
   // Special cases:
   //
   using boost::math::hypergeometric_1F1;
   T tol = boost::math::tools::epsilon<T>() * 200;
   if (std::numeric_limits<T>::digits > std::numeric_limits<double>::digits)
      tol *= 2;
   if (boost::is_class<T>::value)
      tol *= 4;
   // b = 2a
   T computed = hypergeometric_1F1(T(-12.25), T(2 * -12.25), T(6.75));
   T expected = boost::lexical_cast<T>("22.995348157760091167706081204212893687052775606591209203948675272473773725021024450870565197330528784707135828761");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(12.25), T(2 * 12.25), T(6.75));
   expected = boost::lexical_cast<T>("36.47281964229300610642392880149257389834650024065756742702265701321933782423217084029882132197130099355867287657");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-11), T(-12), T(6.75));
   expected = boost::lexical_cast<T>("376.3166426246459656334542608880377435064935064935064935064935064935064935064935064935064935064935064935064935064");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-2), T(-12), T(6.75));
   expected = boost::lexical_cast<T>("2.470170454545454545454545454545454545454545454545454545454545454545454545454545454545454545454545454545454545");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-224), T(-1205), T(6.75));
   expected = boost::lexical_cast<T>("3.497033449657595724636676193024114597507981035316405619832857546161530808157860391434240068189887198094611519953");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(0.5), T(-1205.5), T(-6.75));
   expected = boost::lexical_cast<T>("1.00281149043026925155096279505879868076290060374397866773878698584557482321961231721407215665017657501846692575");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-0.5), T(-1205.5), T(-6.75));
   expected = boost::lexical_cast<T>("0.99719639844965644594352920596780535220516138060108955206195178371227403775248888108818326220977962797312690");
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-12), T(16.25), T(1043.75));
   expected = boost::lexical_cast<T>("1.26527673505477678311707565502355407505496430400394171269315320194708537626079491650410923064978320042481912e20");
   BOOST_CHECK_CLOSE(computed, expected, tol * 3);
   
   computed = hypergeometric_1F1(T(3.5), T(3.5), T(36.25));
   expected = exp(T(36.25));
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-3.5), T(-3.5), T(36.25));
   expected = exp(T(36.25));
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(1), T(2), T(36.25));
   expected = boost::math::expm1(T(36.25)) / T(36.25);
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(10.25), T(9.25), T(36.25));
   expected = exp(T(36.25)) * (T(9.25) + T(36.25)) / T(9.25);
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-10.25), T(-11.25), T(36.25));
   expected = exp(T(36.25)) * (T(-11.25) + T(36.25)) / T(-11.25);
   BOOST_CHECK_CLOSE(computed, expected, tol);
   computed = hypergeometric_1F1(T(-10.25), T(-11.25), T(-36.25));
   expected = exp(T(-36.25)) * (T(-11.25) + T(-36.25)) / T(-11.25);
   BOOST_CHECK_CLOSE(computed, expected, tol);
}

template <class T>
void test_spots(T z, const char* type_name)
{
   test_spots1(z, type_name);
   test_spots2(z, type_name);
   if (std::numeric_limits<T>::digits10 < 20)
      test_spots3(z, type_name);
#ifdef TEST_UNSOLVED
   test_spots4(z, type_name);
#endif
   test_spots5(z, type_name);
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

    Real tol = boost::math::tools::epsilon<Real>() * 5;
    BOOST_CHECK_CLOSE_FRACTION(computed, expected, tol);
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
