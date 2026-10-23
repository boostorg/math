/*
 * Copyright Nick Thompson, 2019
 * Use, modification and distribution are subject to the
 * Boost Software License, Version 1.0. (See accompanying file
 * LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 */

#include "math_unit_test.hpp"
#include <numeric>
#include <utility>
#include <random>
#include <cmath>
#include <limits>
#include <boost/math/tools/test_value.hpp>
#include <boost/math/special_functions/jacobi.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/concepts/real_concept.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

#if __has_include(<stdfloat>)
#  include <stdfloat>
#endif

using std::abs;
using boost::math::jacobi;
using boost::math::jacobi_derivative;

template<typename Real>
void test_to_quadratic()
{
    Real h = 1/Real(8);
    for (Real alpha = -1 + h; alpha < 2; alpha += h) {
        for (Real beta = -1 + h; beta < 2; beta += h) {
            for (Real x = -1; x < 1; x += h) {
                Real expected = 1;
                Real computed = jacobi(0, alpha, beta, x);
                CHECK_ULP_CLOSE(expected, computed, 0);

                expected = (alpha + 1) + (alpha + beta +2)*(x-1)/2;
                computed = jacobi(1, alpha, beta, x);
                CHECK_ULP_CLOSE(expected, computed, 0);

                expected = (alpha + 1)*(alpha+2)/2 + (alpha + 2)*(alpha + beta + 3)*(x-1)/2 + (alpha + beta + 3)*(alpha + beta + 4)*(x-1)*(x-1)/8;
                computed = jacobi(2, alpha, beta, x);
                CHECK_ULP_CLOSE(expected, computed, 1);

            }
        }
    }
}

template<typename Real>
void test_symmetry()
{
    Real h = 1/Real(4);
    for (Real alpha = -1 + h; alpha < 2; alpha += h) {
        for (Real beta = -1 + h; beta < 2; beta += h) {
            for (Real x = -1; x < 1; x += h) {
                for (size_t n = 0; n < 20; n += 2)
                {
                    Real expected = jacobi(n, beta, alpha , -x);
                    Real computed = jacobi(n, alpha, beta, x);
                    CHECK_ULP_CLOSE(expected, computed, 0);

                    expected = jacobi(n+1, beta, alpha, -x);
                    computed = -jacobi(n+1, alpha, beta, x);
                    CHECK_ULP_CLOSE(expected, computed, 0);
                }
            }
        }
    }
}

template<typename Real>
void test_derivative()
{
    Real h = 1/Real(4);
    for (Real alpha = -1 + h; alpha < 2; alpha += h) {
        for (Real beta = -1 + h; beta < 2; beta += h) {
            for (Real x = -1; x < 1; x += h) {
                Real expected = 0;
                Real computed = jacobi_derivative(0, alpha, beta, x, 1);
                CHECK_ULP_CLOSE(expected, computed, 0);

                expected = (alpha + beta + 2)/2;
                computed = jacobi_derivative(1, alpha, beta, x, 1);
                CHECK_ULP_CLOSE(expected, computed, 0);

                expected = (alpha + 2)*(alpha + beta + 3)/2 + (alpha + beta + 3)*(alpha + beta + 4)*(x-1)/4;
                computed = jacobi_derivative(2, alpha, beta, x, 1);
                CHECK_ULP_CLOSE(expected, computed, 0);

                expected = (alpha + beta + 3)*(alpha + beta + 4)/4;
                computed = jacobi_derivative(2, alpha, beta, x, 2);
                CHECK_ULP_CLOSE(expected, computed, 0);
            }
        }
    }
}

template <class T>
void test_zeros(T tol, const char* t)
{
   std::cout << "Testing zeros of Jacobi polynomials for type " << t << std::endl;

   unsigned N = 10;
   T alpha = -static_cast<T>(0.5);
   T beta = -static_cast<T>(0.25);
   std::vector<T> roots = boost::math::jacobi_zeros<T>(N, alpha, beta);

   // Values from JacobiRoots.txt, truncated to the original number of digits.
   CHECK_ABSOLUTE_ERROR(roots[0], BOOST_MATH_TEST_VALUE(T, -0.98043481087135728948225748512377413191431377879335791867734106797588116008793062549212633665088501425444157046969), tol);
   CHECK_ABSOLUTE_ERROR(roots[1], BOOST_MATH_TEST_VALUE(T, -0.87471795355816395604640367310993598116332684497631361442801627612702557667163210395209101070585349427142111164288), tol);
   CHECK_ABSOLUTE_ERROR(roots[2], BOOST_MATH_TEST_VALUE(T, -0.68547729631419309000494293417916787382452464017484521300870391546908609386112534316875175649816490820761997344333), tol);
   CHECK_ABSOLUTE_ERROR(roots[3], BOOST_MATH_TEST_VALUE(T, -0.43073485032125228809815001045027898152469157970218561300262817883164860859175888263649327736176401127435226007673), tol);
   CHECK_ABSOLUTE_ERROR(roots[4], BOOST_MATH_TEST_VALUE(T, -0.13481395267014544369612252045627773208063462409759759746513277127211221648854974876657289823245618174333534107812), tol);
   CHECK_ABSOLUTE_ERROR(roots[5], BOOST_MATH_TEST_VALUE(T, 0.174025036041258702833387110360958607634851774908859934577740491653957330826514516755665614613877586968329344222486), tol);
   CHECK_ABSOLUTE_ERROR(roots[6], BOOST_MATH_TEST_VALUE(T, 0.466287142403906575515711325793441163303393780413881146080881185495212690974675349190272063596260412183693498775484), tol);
   CHECK_ABSOLUTE_ERROR(roots[7], BOOST_MATH_TEST_VALUE(T, 0.714060291625917422384478107667342985650985129021464154416551477200455821467118506327105981403288530428645575696712), tol);
   CHECK_ABSOLUTE_ERROR(roots[8], BOOST_MATH_TEST_VALUE(T, 0.893681191657027897970840390419389902657653328912723300666540830669739848619590156972674550185798367491416410935830), tol);
   CHECK_ABSOLUTE_ERROR(roots[9], BOOST_MATH_TEST_VALUE(T, 0.987995331877131338753329818948431911390477584357501290969978354526517833942968304640446939779768842549215297210134), tol);
}

template <class T>
void test_zero_special_case()
{
   std::cout << "Testing special cases for zeros" << std::endl;

   CHECK_THROW(boost::math::jacobi_zeros<T>(0, 1, 1), std::domain_error);
   CHECK_THROW(boost::math::jacobi_zeros<T>(10, -1, 1), std::domain_error);
   CHECK_THROW(boost::math::jacobi_zeros<T>(10, 1, -1), std::domain_error);
   CHECK_THROW(boost::math::jacobi_zeros<T>(10, std::numeric_limits<T>::quiet_NaN(), 1), std::domain_error);
   CHECK_THROW(boost::math::jacobi_zeros<T>(10, 1, std::numeric_limits<T>::quiet_NaN()), std::domain_error);
   CHECK_THROW(boost::math::jacobi_zeros<T>(std::numeric_limits<unsigned>::quiet_NaN(), 1, 1), std::domain_error);

   using ignore_policy = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error>,
                                                       boost::math::policies::evaluation_error<boost::math::policies::ignore_error> >;
   CHECK_NAN((boost::math::jacobi_zeros<T, ignore_policy>(0, 1, 1)[0]));
   CHECK_NAN((boost::math::jacobi_zeros<T, ignore_policy>(std::numeric_limits<unsigned>::quiet_NaN(), 1, 1)[0]));

   std::vector<T> zeros = boost::math::jacobi_zeros<T>(1, 0, 0);
   CHECK_EQUAL(zeros.size(), static_cast<std::size_t>(1));
   CHECK_EQUAL(zeros[0], static_cast<T>(0));
}

int main()
{
    #ifdef __STDCPP_FLOAT32_T__
    
    test_symmetry<std::float32_t>();
    test_derivative<std::float32_t>();

    #else

    test_symmetry<float>();
    test_derivative<float>();

    #endif
    
    #ifdef __STDCPP_FLOAT64_T__

    test_to_quadratic<std::float64_t>();
    test_symmetry<std::float64_t>();
    test_derivative<std::float64_t>();

    #else

    test_to_quadratic<double>();
    test_symmetry<double>();
    test_derivative<double>();

    #endif

    #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
    
    test_to_quadratic<long double>();
    test_symmetry<long double>();
    test_derivative<long double>();
    
    #endif

    #ifdef BOOST_HAS_FLOAT128
    
    test_to_quadratic<boost::multiprecision::float128>();
    test_symmetry<boost::multiprecision::float128>();
    test_derivative<boost::multiprecision::float128>();
    
    #endif

    test_zeros(2*std::numeric_limits<float>::epsilon(), "float");
    test_zeros(std::numeric_limits<double>::epsilon(), "double");
    test_zeros(std::numeric_limits<long double>::epsilon(), "long double");

    using cpp_bin_float_100 = boost::multiprecision::cpp_bin_float_100;
    test_zeros<cpp_bin_float_100>(std::numeric_limits<cpp_bin_float_100>::epsilon(), "cpp_bin_float_100");
    #ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
    test_zeros(boost::math::concepts::real_concept(1e-15), "real_concept");
    #endif

    test_zero_special_case<double>();

    return boost::math::test::report_errors();
}
