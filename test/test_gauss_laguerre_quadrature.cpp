// Copyright Jacob Hass, 2026
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
      Real Q = gauss_laguerre<Real, N>().integrate(even, nullptr, &L1);
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
      CHECK_ULP_CLOSE(static_cast<Real>(table::abscissa()[i]), on_demand::abscissa()[i], 10);
      CHECK_ULP_CLOSE(static_cast<Real>(table::weights()[i]), on_demand::weights()[i], 32);
   }
}

// Integral of exp(ix) exp(-x) over the real line is (0.5, 0.5).
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
   using vector = std::valarray<double>;
   const vector zero(0.0, 3);
   auto norm = [](const vector& v) { return std::sqrt((v * v).sum()); };
   auto f = [](double x) { return vector{1.0, x * x, std::cos(x)}; };
   double L1;
   vector Q = gauss_laguerre<double, N>().integrate(f, zero, norm, nullptr, &L1);
   CHECK_ULP_CLOSE(1.0, Q[0], 10);
   CHECK_ULP_CLOSE(2.0, Q[1], 20);
   CHECK_ULP_CLOSE(0.5, Q[2], 16);
   CHECK_LE(norm(Q), L1);
   // The estimate is the norm of the difference from the anti-Gauss rule, at least that of any component.
   double error, L1_error, cos_error;
   const vector Q_error = gauss_laguerre<double, N>().integrate(f, zero, norm, &error, &L1_error);
   CHECK_ULP_CLOSE(Q_error[2], Q[2], 2);
   CHECK_ULP_CLOSE(L1_error, L1, 2);
   gauss_laguerre<double, N>().integrate([](double x) { return std::cos(x); }, &cos_error, nullptr);
   CHECK_LE(cos_error, error);
   CHECK_LE(2 * std::numeric_limits<double>::epsilon() * L1, error);
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

// Past N = 706 in double the recurrence overflows, which must be reported through the policy.
void test_overflow()
{
   using ignore_policy = boost::math::policies::policy<boost::math::policies::evaluation_error<boost::math::policies::ignore_error> >;
   double L1 = 0;
   CHECK_NAN((gauss_laguerre<double, 1000, ignore_policy>().integrate([](double x) { return std::cos(x); }, nullptr, &L1)));
   CHECK_NAN(L1);
   using vector = std::valarray<double>;
   vector Q = gauss_laguerre<double, 1000, ignore_policy>().integrate([](double x) { return vector{1.0, x}; }, vector(0.0, 2),
                                                                    [](const vector& v) { return std::sqrt((v * v).sum()); });
   CHECK_EQUAL(Q.size(), std::size_t(2));
   CHECK_NAN(Q[0]);
   CHECK_NAN(Q[1]);
   double error = 0;
   gauss_laguerre<double, 1000, ignore_policy>().integrate([](double x) { return std::cos(x); }, &error, nullptr);
   CHECK_NAN(error);
#ifndef BOOST_NO_EXCEPTIONS
   CHECK_THROW((gauss_laguerre<double, 1000>().integrate([](double x) { return std::cos(x); })), boost::math::evaluation_error);
#endif
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

// For smooth f the estimate is close to the true error until that reaches rounding.
template <class Real, unsigned N>
void test_error_estimate()
{
   using std::abs;
   using std::cos;
   using std::sqrt;
   using rule = gauss_laguerre<Real, N>;
   auto f = [](const Real& x) { return cos(x); };
   Real error, L1, L1_alone;
   const Real Q = rule().integrate(f, &error, &L1);
   CHECK_ULP_CLOSE(Q, rule().integrate(f, nullptr, &L1_alone), 2);
   CHECK_ULP_CLOSE(L1, L1_alone, 2);
   CHECK_ULP_CLOSE(Q, rule().integrate(f, nullptr, nullptr), 2);
   const Real true_error = abs(Q - Real(0.5));
   if (true_error > 1000 * std::numeric_limits<Real>::epsilon() * L1)
   {
      CHECK_LE(error, 2 * true_error);
      CHECK_LE(true_error, 2 * error);
   }
   CHECK_LE(2 * std::numeric_limits<Real>::epsilon() * L1, error);
   // The integral of exp(-x) / (1 + x) is the Gompertz constant; the pole at -1 makes the error fall slowly.
   const Real Qg = rule().integrate([](const Real& x) { return 1 / (1 + x); }, &error, nullptr);
   const Real true_error_g = abs(Qg - static_cast<Real>(0.5963473623231940743410784993692793760742L));
   CHECK_LE(error, 2 * true_error_g);
   CHECK_LE(true_error_g, 2 * error);
   const Real Qh = rule().integrate([](const Real& x) { return sqrt(1 + x); }, &error, nullptr);
   const Real true_error_h = abs(Qh - static_cast<Real>(1.378936078070656053021675619957108958174L));
   CHECK_LE(error, 2 * true_error_h);
   CHECK_LE(true_error_h, 2 * error);
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
      const Real Q = rule.integrate([k](const Real& x) { return pow(x, Real(k)); }, nullptr, &L1);
      CHECK_ULP_CLOSE(boost::math::tgamma(k + alpha + 1), Q, static_cast<int>(moment_tol * (k + 1)));
   }
   Real error, L1;
   const Real Q = rule.integrate([](const Real& x) { return cos(x); }, &error, &L1);
   const Real expected = boost::math::tgamma(alpha + 1) * cos((alpha + 1) * boost::math::constants::pi<Real>() / 4) / pow(Real(2), (alpha + 1) / 2);
   const Real true_error = abs(Q - expected);
   // For fewer points the estimate is preasymptotic: 8 times the error for N = 2 and alpha = 2.5.
   if ((N >= 12) && (true_error > 1000 * std::numeric_limits<Real>::epsilon() * L1))
   {
      CHECK_LE(error, 2 * true_error);
      CHECK_LE(true_error, 2 * error);
   }
   CHECK_LE(2 * std::numeric_limits<Real>::epsilon() * L1, error);
}

void test_alpha_domain()
{
   using ignore_policy = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error> >;
   for (double alpha : { -1.0, -2.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity() })
   {
      double error = 0;
      CHECK_NAN((gauss_laguerre<double, 5, ignore_policy>(alpha).integrate([](double x) { return x; }, &error, nullptr)));
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
   // For alpha near -1 the smallest abscissa, and so the largest weight, is less accurate.
   for (double alpha : { -0.9, -0.5, 0.5, 2.5 })
   {
      const double tol = alpha < 0 ? 64 : 16;
      test_alpha<double, 2>(alpha, tol);
      test_alpha<double, 5>(alpha, tol);
      test_alpha<double, 12>(alpha, tol);
      test_alpha<double, 20>(alpha, tol);
      test_alpha<double, 30>(alpha, tol);
   }
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_alpha<long double, 12>(-0.5L, 64);
#endif
   test_alpha_domain();

   return boost::math::test::report_errors();
}
