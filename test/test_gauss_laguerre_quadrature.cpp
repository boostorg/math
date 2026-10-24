// Copyright Jacob Hass, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include "gauss_quadrature_test_helpers.hpp"
#include <cmath>
#include <complex>
#include <limits>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/gauss_laguerre.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/special_functions/factorials.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

using boost::math::quadrature::gauss_laguerre;
using boost::math::quadrature::detail::gauss_laguerre_detail;

// An N-point rule integrates x^(n) exp(-x) exactly for n <= 2N - 1, giving Gamma(n + 1).
// The tolerance grows with n because x^(n) magnifies the rounding of the largest abscissa
// by a factor of n.
template <class Real, unsigned N>
void test_moments()
{
   // Gamma(k + 1) by its recurrence
   Real expected;
   for (unsigned k = 0; k < 2 * N - 1; ++k)
   {
      std::cout << "Testing N = " << N << ", k = " << k << std::endl;
      expected = (k == 0) ? Real(1) : boost::math::factorial<Real>(k);
      auto even = [k](const Real& x) { Real r = 1; for (unsigned j = 0; j < k; ++j) r *= x; return r; };
      Real L1;
      Real Q = gauss_laguerre<Real, N>().integrate(even, &L1);
      std::cout << std::setprecision(17) << expected << " " << Q << " " << L1 << std::endl;
      CHECK_ULP_CLOSE(expected, Q, 16 * (k + 1));
      // The integrand is nonnegative, so the L1 norm is the integral itself.
      CHECK_ULP_CLOSE(Q, L1, 16 * (k + 1));
   }
}

// The tabulated abscissas and weights must agree with those computed on demand.
template <class Real, unsigned N, unsigned Category>
void test_tables_match_on_demand()
{
   using table = gauss_laguerre_detail<Real, N, Category>;
   using on_demand = gauss_laguerre_detail<Real, N, 999>;
   CHECK_EQUAL(table::abscissa().size(), on_demand::abscissa().size());
   for (unsigned i = 0; i < table::abscissa().size(); ++i)
   {
      CHECK_ULP_CLOSE(static_cast<Real>(table::abscissa()[i]), on_demand::abscissa()[i], 2);
      // cpp_bin_float's exp, which the weights computed on demand take at the nodes, is out by 4 ulps at x = 19.
      CHECK_ULP_CLOSE(static_cast<Real>(table::weights()[i]), on_demand::weights()[i], 16);
   }
}

// Integral of exp(ix) exp(-x) over [0, infinity) is (0.5, 0.5).
template <class Complex, unsigned N>
void test_complex(typename Complex::value_type tol)
{
   using Real = typename Complex::value_type;
   using std::cos;
   using std::exp;
   using std::sin;
   auto f = [](const Real& x) { return Complex(cos(x), sin(x)); };
   Complex Q = gauss_laguerre<Real, N>().integrate(f);
   CHECK_ABSOLUTE_ERROR(Real(0.5), Q.real(), tol);
   CHECK_ABSOLUTE_ERROR(Real(0.5), Q.imag(), tol);
}

template <unsigned N>
void test_vector_valued()
{
   check_vector_valued(gauss_laguerre<double, N>(), 1.0, 2.0, 0.5, 10, 20);
}

// Evaluating L_N directly overflows long before these N; the orthonormal recurrence does not.
template <unsigned N>
void test_large_N()
{
   using rule = gauss_laguerre<double, N>;
   CHECK_EQUAL(rule().abscissa().size(), static_cast<std::size_t>(N));
   for (unsigned i = 1; i < rule().abscissa().size(); ++i)
      CHECK_LE(rule().abscissa()[i - 1], rule().abscissa()[i]);
   CHECK_ULP_CLOSE(1.0, rule().integrate([](double) { return 1.0; }), 52);
   CHECK_ULP_CLOSE(0.5, rule().integrate([](double x) { return std::cos(x); }), 52);
}

// Past N = 326 the recurrence overflows in double. The Gauss rule comes from the differential equation, which does not,
// Past N = 326 the recurrence overflows in double, but the differential equation does not: neither rule overflows.
void test_overflow()
{
   double L1 = 0, error = 0;
   CHECK_ULP_CLOSE(0.5, (gauss_laguerre<double, 1000>().integrate([](double x) { return std::cos(x); }, &L1)), 64);
   CHECK_ULP_CLOSE(0.5, (gauss_laguerre<double, 1000>().integrate_averaged([](double x) { return std::cos(x); }, &error, &L1)), 64);
   CHECK_LE(error, 64 * std::numeric_limits<double>::epsilon());
}

// The anti-Gauss rule is the Gauss rule of 2I - G_N, so on polynomials of degree 2N + 1 its error is that of G_N
// with the opposite sign. Its N + 1 abscissas interlace those of G_N and are positive (Laurie 1996, Theorem 2),
// and its weights are positive and sum to Gamma(alpha + 1).
template <class Real>
void test_anti_gauss(unsigned N, Real alpha = 0)
{
   using std::abs;
   using family = boost::math::detail::laguerre_family<Real>;
   const boost::math::detail::orthogonal_polynomial<Real, family> gauss(N, family(N, alpha));
   const boost::math::detail::orthogonal_polynomial<Real, family> anti(N + 1, family(N + 1, alpha), true);
   const std::vector<Real>& x = anti.abscissa();
   const std::vector<Real>& w = anti.weights();
   CHECK_EQUAL(x.size(), static_cast<std::size_t>(N + 1));
   CHECK_LE(std::numeric_limits<Real>::min(), x[0]);
   Real sum = 0;
   for (std::size_t i = 0; i < x.size(); ++i)
   {
      CHECK_LE(Real(0), w[i]);
      sum += w[i];
   }
   CHECK_ULP_CLOSE(boost::math::tgamma(alpha + 1), sum, 8);
   for (std::size_t i = 0; i < N; ++i)
   {
      CHECK_LE(x[i], gauss.abscissa()[i]);
      CHECK_LE(gauss.abscissa()[i], x[i + 1]);
   }
   // p = sum c_k L_k^alpha with arbitrary c_k, k <= 2N + 1, whose integral is c_0 Gamma(alpha + 1). The error of G_N
   // on L_2N is 1 / binomial(2N, N), so this only constrains the rule for small N; rounding grows with N.
   std::vector<Real> c(2 * N + 2);
   for (unsigned k = 0; k < c.size(); ++k)
   {
      c[k] = Real((k * 7919) % 17) / 8 - 1;
   }
   const Real exact = c[0] * boost::math::tgamma(alpha + 1);
   auto p = [&c, alpha](const Real& t)
   {
      Real r = c[0], L = 1, L_previous = 0;
      for (unsigned k = 1; k < c.size(); ++k)
      {
         Real L_next = ((2 * k - 1 + alpha - t) * L - (k - 1 + alpha) * L_previous) / k;
         L_previous = L;
         L = L_next;
         r += c[k] * L;
      }
      return r;
   };
   Real gauss_result = 0, anti_result = 0, L1 = 0;
   for (std::size_t i = 0; i < N; ++i)
   {
      gauss_result += gauss.weights()[i] * p(gauss.abscissa()[i]);
      L1 += gauss.weights()[i] * abs(p(gauss.abscissa()[i]));
   }
   for (std::size_t i = 0; i < x.size(); ++i)
   {
      anti_result += w[i] * p(x[i]);
      L1 += w[i] * abs(p(x[i]));
   }
   CHECK_ABSOLUTE_ERROR(Real(0), (exact - gauss_result) + (exact - anti_result), (N + 8) * std::numeric_limits<Real>::epsilon() * L1);
}

template <class Real, unsigned N>
void test_error_estimate()
{
   using std::cos;
   using std::sqrt;
   check_error_estimate(gauss_laguerre<Real, N>(), [](const Real& x) { return cos(x); }, Real(0.5));
   // The integral of exp(-x) / (1 + x) is the Gompertz constant; the pole at -1 makes the error fall slowly.
   check_error_estimate(gauss_laguerre<Real, N>(), [](const Real& x) { return 1 / (1 + x); }, static_cast<Real>(0.5963473623231940743410784993692793760742L));
   check_error_estimate(gauss_laguerre<Real, N>(), [](const Real& x) { return sqrt(1 + x); }, static_cast<Real>(1.378936078070656053021675619957108958174L));
}

// The rule for x^alpha exp(-x) integrates x^k exactly for k <= 2N - 1, giving Gamma(k + alpha + 1), and its error
// estimate tracks the true error for the integral of x^alpha exp(-x) cos(x), Gamma(alpha + 1) cos((alpha + 1) pi / 4) / 2^((alpha + 1) / 2).
template <class Real, unsigned N>
void test_alpha(Real alpha, Real moment_tol)
{
   using std::abs;
   using std::cos;
   using std::pow;
   const gauss_laguerre<Real, N> rule(alpha);
   for (unsigned k = 0; k < 2 * N; ++k)
   {
      Real L1;
      const Real Q = rule.integrate([k](const Real& x) { return pow(x, Real(k)); }, &L1);
      CHECK_ULP_CLOSE(boost::math::tgamma(k + alpha + 1), Q, static_cast<int>(moment_tol * (k + 1)));
   }
   Real error, L1;
   const Real Q = rule.integrate_averaged([](const Real& x) { return cos(x); }, &error, &L1);
   const Real expected = boost::math::tgamma(alpha + 1) * cos((alpha + 1) * boost::math::constants::pi<Real>() / 4) / pow(Real(2), (alpha + 1) / 2);
   const Real gauss_error = abs(rule.integrate([](const Real& x) { return cos(x); }) - expected);
   // For fewer points the estimate is preasymptotic: 8 times the error for N = 2 and alpha = 2.5, and G and A
   // need not bracket the integral.
   if ((N >= 12) && (gauss_error > 1000 * std::numeric_limits<Real>::epsilon() * L1))
   {
      CHECK_LE(abs(Q - expected), error);
      CHECK_LE(error, 2 * gauss_error);
      CHECK_LE(gauss_error, 2 * error);
   }
   CHECK_LE(2 * std::numeric_limits<Real>::epsilon() * L1, error);
}

template <unsigned N>
void test_anti_gauss_accessors(double alpha)
{
   using family = boost::math::detail::laguerre_family<double>;
   const gauss_laguerre<double, N> rule(alpha);
   const boost::math::detail::orthogonal_polynomial<double, family> anti(N + 1, family(N + 1, alpha), true);
   CHECK_EQUAL(rule.anti_gauss_abscissa().size(), static_cast<std::size_t>(N + 1));
   for (unsigned i = 0; i <= N; ++i)
   {
      CHECK_EQUAL(anti.abscissa()[i], rule.anti_gauss_abscissa()[i]);
      CHECK_EQUAL(anti.weights()[i], rule.anti_gauss_weights()[i]);
   }
   // The averaged rule is exact on x^(2N + 1), whose integral is Gamma(2N + 2 + alpha).
   auto f = [](double x) { return std::pow(x, 2 * N + 1); };
   CHECK_ULP_CLOSE(boost::math::tgamma(2 * N + 2 + alpha), rule.integrate_averaged(f), 64);
}

void test_alpha_domain()
{
   using ignore_policy = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error> >;
   for (double alpha : { -1.0, -2.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() })
   {
      double error = 0;
      CHECK_NAN((gauss_laguerre<double, 5, ignore_policy>(alpha).integrate([](double x) { return x; })));
      CHECK_NAN((gauss_laguerre<double, 5, ignore_policy>(alpha).integrate_averaged([](double x) { return x; }, &error)));
      CHECK_NAN(error);
#ifndef BOOST_NO_EXCEPTIONS
      CHECK_THROW((gauss_laguerre<double, 5>(alpha)), std::domain_error);
#endif
   }
}

template <class Real>
void test_moments_all_N()
{
   test_moments<Real, 2>();
   test_moments<Real, 5>();
   test_moments<Real, 7>();
   test_moments<Real, 10>();
   test_moments<Real, 15>();
   test_moments<Real, 20>();
}

int main()
{
   // Factorial of floats overflow for larger N > 10
   test_moments<float, 2>();
   test_moments<float, 5>();
   test_moments<float, 7>();
   test_moments<float, 10>();

   test_moments_all_N<double>();
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_moments_all_N<long double>();
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_moments<boost::math::concepts::real_concept, 7>();
   test_moments<boost::math::concepts::real_concept, 10>();
#endif
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

   test_complex<std::complex<double>, 30>(50 * std::numeric_limits<double>::epsilon());
   test_complex<boost::multiprecision::cpp_complex_quad, 60>(50 * std::numeric_limits<boost::multiprecision::cpp_complex_quad::value_type>::epsilon());

   test_vector_valued<30>();
   test_vector_valued<40>();

   // Overflow occurs at about N = 200
   test_large_N<100>();
   test_large_N<190>();

   test_overflow();

   for (unsigned n = 1; n <= 50; ++n)
   {
      test_anti_gauss<double>(n);
      test_anti_gauss<double>(n, -0.9);
      test_anti_gauss<double>(n, -0.99);
      test_anti_gauss<double>(n, 2.5);
   }
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_anti_gauss<long double>(7);
   test_anti_gauss<long double>(20);
#endif
   test_error_estimate<double, 2>();
   test_error_estimate<double, 5>();
   test_error_estimate<double, 10>();
   test_error_estimate<double, 20>();
   test_error_estimate<double, 40>();
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_error_estimate<long double, 15>();
#endif
   for (double alpha : { -0.9, -0.5, 0.5, 2.5 })
   {
      const double tol = 4;
      test_alpha<double, 2>(alpha, tol);
      test_alpha<double, 5>(alpha, tol);
      test_alpha<double, 12>(alpha, tol);
      test_alpha<double, 20>(alpha, tol);
      test_alpha<double, 30>(alpha, tol);
   }
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_alpha<long double, 12>(-0.5L, 4);
#endif
   test_alpha_domain();
   test_anti_gauss_accessors<7>(0);
   test_anti_gauss_accessors<10>(0.5);

   return boost::math::test::report_errors();
}
