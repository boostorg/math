// Copyright Jacob Hass, 2026
// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <valarray>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/gauss_hermite.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

using boost::math::quadrature::gauss_hermite;
using boost::math::quadrature::detail::gauss_hermite_detail;

// An N-point rule integrates x^(2k) exp(-x^2) exactly for 2k <= 2N - 1, giving Gamma(k + 1/2).
// Odd moments vanish by symmetry. The tolerance grows with k because x^(2k) magnifies the
// rounding of the largest abscissa by a factor of 2k.
template <class Real, unsigned N>
void test_moments()
{
   // Gamma(k + 1/2) by its recurrence from sqrt(pi), which is more accurate than tgamma for some types.
   Real expected = boost::math::constants::root_pi<Real>();
   for (unsigned k = 0; k < N; ++k)
   {
      if (k > 0)
         expected *= Real(k) - Real(0.5);
      auto even = [k](const Real& x) { Real r = 1; for (unsigned j = 0; j < 2 * k; ++j) r *= x; return r; };
      Real L1;
      Real Q = gauss_hermite<Real, N>::integrate(even, &L1);
      CHECK_ULP_CLOSE(expected, Q, 16 * (k + 1));
      // The integrand is nonnegative, so the L1 norm is the integral itself.
      CHECK_ULP_CLOSE(Q, L1, 16 * (k + 1));

      auto odd = [k](const Real& x) { Real r = x; for (unsigned j = 0; j < 2 * k; ++j) r *= x; return r; };
      CHECK_EQUAL(Real(gauss_hermite<Real, N>::integrate(odd)), Real(0));
   }
}

// The tabulated abscissas and weights must agree with those computed on demand.
template <class Real, unsigned N, unsigned Category>
void test_tables_match_on_demand()
{
   using table = gauss_hermite_detail<Real, N, Category>;
   using on_demand = gauss_hermite_detail<Real, N, 999>;
   CHECK_EQUAL(table::abscissa().size(), on_demand::abscissa().size());
   for (unsigned i = 0; i < table::abscissa().size(); ++i)
   {
      CHECK_ULP_CLOSE(static_cast<Real>(table::abscissa()[i]), on_demand::abscissa()[i], 4);
      CHECK_ULP_CLOSE(static_cast<Real>(table::weights()[i]), on_demand::weights()[i], 32);
   }
}

// Integral of exp(ix) exp(-x^2) over the real line is sqrt(pi) exp(-1/4).
template <class Complex, unsigned N>
void test_complex(std::size_t ulps, typename Complex::value_type imaginary_tol)
{
   using Real = typename Complex::value_type;
   using std::cos;
   using std::exp;
   using std::sin;
   auto f = [](const Real& x) { return Complex(cos(x), sin(x)); };
   Complex Q = gauss_hermite<Real, N>::integrate(f);
   Real expected = boost::math::constants::root_pi<Real>() * exp(-Real(1) / 4);
   CHECK_ULP_CLOSE(expected, Q.real(), ulps);
   CHECK_ABSOLUTE_ERROR(Real(0), Q.imag(), imaginary_tol);
}

template <unsigned N>
void test_vector_valued()
{
   using vector = std::valarray<double>;
   const vector zero(0.0, 3);
   auto norm = [](const vector& v) { return std::sqrt((v * v).sum()); };
   auto f = [](double x) { return vector{1.0, x * x, std::cos(x)}; };
   double L1;
   vector Q = gauss_hermite<double, N>::integrate(f, zero, norm, &L1);
   const double root_pi = boost::math::constants::root_pi<double>();
   CHECK_ULP_CLOSE(root_pi, Q[0], 4);
   CHECK_ULP_CLOSE(root_pi / 2, Q[1], 4);
   CHECK_ULP_CLOSE(root_pi * std::exp(-0.25), Q[2], 16);
   CHECK_LE(norm(Q), L1);
}

// Evaluating H_N directly overflows long before these N; the orthonormal recurrence does not.
template <unsigned N>
void test_large_N()
{
   using rule = gauss_hermite<double, N>;
   const double root_pi = boost::math::constants::root_pi<double>();
   CHECK_EQUAL(rule::abscissa().size(), static_cast<std::size_t>((N + 1) / 2));
   for (unsigned i = 1; i < rule::abscissa().size(); ++i)
      CHECK_LE(rule::abscissa()[i - 1], rule::abscissa()[i]);
   CHECK_ULP_CLOSE(root_pi, rule::integrate([](double) { return 1.0; }), 16);
   CHECK_ULP_CLOSE(root_pi * std::exp(-0.25), rule::integrate([](double x) { return std::cos(x); }), 16);
}

// Past N = 706 in double the recurrence overflows, which must be reported through the policy.
void test_overflow()
{
   using ignore_policy = boost::math::policies::policy<boost::math::policies::evaluation_error<boost::math::policies::ignore_error> >;
   double L1 = 0;
   CHECK_NAN((gauss_hermite<double, 1000, ignore_policy>::integrate([](double x) { return std::cos(x); }, &L1)));
   CHECK_NAN(L1);
   using vector = std::valarray<double>;
   vector Q = gauss_hermite<double, 1000, ignore_policy>::integrate([](double x) { return vector{1.0, x}; }, vector(0.0, 2),
                                                                    [](const vector& v) { return std::sqrt((v * v).sum()); });
   CHECK_EQUAL(Q.size(), std::size_t(2));
   CHECK_NAN(Q[0]);
   CHECK_NAN(Q[1]);
#ifndef BOOST_NO_EXCEPTIONS
   CHECK_THROW((gauss_hermite<double, 1000>::integrate([](double x) { return std::cos(x); })), boost::math::evaluation_error);
#endif
}

template <class Real>
void test_moments_all_N()
{
   test_moments<Real, 1>();
   test_moments<Real, 2>();
   test_moments<Real, 5>();
   test_moments<Real, 7>();
   test_moments<Real, 10>();
   test_moments<Real, 15>();
   test_moments<Real, 20>();
}

int main()
{
   test_moments_all_N<float>();
   test_moments_all_N<double>();
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_moments_all_N<long double>();
#endif
   test_moments<boost::math::concepts::real_concept, 7>();
   test_moments<boost::math::concepts::real_concept, 10>();
   test_moments_all_N<boost::multiprecision::cpp_bin_float_50>();
   test_moments<boost::multiprecision::cpp_bin_float_quad, 15>();
   test_moments<boost::multiprecision::cpp_bin_float_quad, 20>();
#ifdef BOOST_HAS_FLOAT128
   test_moments<boost::multiprecision::float128, 15>();
   test_moments<boost::multiprecision::float128, 20>();
#endif

   test_tables_match_on_demand<double, 7, 0>();
   test_tables_match_on_demand<double, 10, 0>();
   test_tables_match_on_demand<double, 15, 0>();
   using cpp_bin_float_100 = boost::multiprecision::cpp_bin_float_100;
   test_tables_match_on_demand<cpp_bin_float_100, 7, 4>();
   test_tables_match_on_demand<cpp_bin_float_100, 10, 4>();
   test_tables_match_on_demand<cpp_bin_float_100, 15, 4>();

   test_complex<std::complex<double>, 15>(4, 4 * std::numeric_limits<double>::epsilon());
   test_complex<boost::multiprecision::cpp_complex_quad, 30>(64, 1e-32);

   test_vector_valued<15>();
   test_vector_valued<20>();

   test_large_N<200>();
   test_large_N<500>();
   test_large_N<706>();
   test_overflow();

   return boost::math::test::report_errors();
}
