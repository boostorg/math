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
#include <boost/math/tools/test_value.hpp>
#include <boost/math/concepts/real_concept.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/gauss_jacobi.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/special_functions/factorials.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

using boost::math::quadrature::gauss_jacobi;

// Gauss-Jacobi integration is exact for polynomials of degree less than 2N - 1
// We use a polynomial of degree 10, so should be exact for N >= 6
template <class Real, unsigned N>
void test_moments(Real tol, const char* type_name)
{
   using std::pow;
   std::cout << "Testing Moments for N = " << N << " and type " << type_name << std::endl;
   auto f = [](const Real& x) { return pow(x, Real(10.0)); };
   Real L1;
   gauss_jacobi<Real> integrator = gauss_jacobi<Real>(N, -Real(0.5), -Real(0.25));
   Real Q = integrator.integrate(f, nullptr, &L1);

   CHECK_ABSOLUTE_ERROR(Q, Real(0.604206393406616872500520394161000920557127865245), tol);
   // The integrand is nonnegative, so the L1 norm is the integral itself.
   CHECK_ABSOLUTE_ERROR(Q, L1, tol);
}

template <class Complex, unsigned N>
void test_complex(typename Complex::value_type tol)
{
   using Real = typename Complex::value_type;
   using std::cos;
   using std::exp;
   using std::sin;
   auto f = [](const Real& x) { return Complex(cos(x), sin(x)); };
   gauss_jacobi<Real> integrator = gauss_jacobi<Real>(N, 1.0, 1.0);
   Complex Q = integrator.integrate(f);
   Real real_expected = Real(4) * sin(Real(1)) - Real(4) * cos(Real(1));
   CHECK_ABSOLUTE_ERROR(real_expected, Q.real(), tol);
   CHECK_ABSOLUTE_ERROR(Real(0.0), Q.imag(), tol);
}

template <unsigned N>
void test_vector_valued()
{
   using vector = std::valarray<double>;
   const vector zero(0.0, 3);
   auto norm = [](const vector& v) { return std::sqrt((v * v).sum()); };
   auto f = [](double x) { return vector{1.0, x * x, std::cos(x)}; };
   double L1;
   gauss_jacobi<double> integrator = gauss_jacobi<double>(N, 1.0, 1.0);
   vector Q = integrator.integrate(f, zero, norm, nullptr, &L1);
   CHECK_ULP_CLOSE(4.0 / 3.0, Q[0], 10);
   CHECK_ULP_CLOSE(4.0 / 15.0, Q[1], 20);
   CHECK_ULP_CLOSE(4.0 * sin(1.0) - 4.0 * cos(1.0), Q[2], 16);
   CHECK_LE(norm(Q), L1);
   // The estimate is the norm of the difference from the anti-Gauss rule, at least that of any component.
   double error, L1_error, cos_error;
   const vector Q_error = integrator.integrate(f, zero, norm, &error, &L1_error);
   CHECK_ULP_CLOSE(Q_error[2], Q[2], 2);
   CHECK_ULP_CLOSE(L1_error, L1, 2);
   integrator.integrate([](double x) { return std::cos(x); }, &cos_error, nullptr);
   CHECK_LE(cos_error, error);
   CHECK_LE(2 * std::numeric_limits<double>::epsilon() * L1, error);
}

// The anti-Gauss rule is the Gauss rule of 2I - G_N, so on polynomials of degree 2N + 1 its error is that of G_N
// with the opposite sign; its N + 1 abscissas interlace those of G_N, and its weights are positive and sum to
// 2^(alpha + beta + 1) B(alpha + 1, beta + 1).
template <class Real>
void test_anti_gauss(unsigned N, Real alpha, Real beta)
{
   using std::abs;
   using std::pow;
   using family = boost::math::detail::jacobi_family<Real>;
   gauss_jacobi<Real> rule(N, alpha, beta);
   const boost::math::detail::orthogonal_polynomial<Real, family> anti(N + 1, family(N + 1, alpha, beta), true);
   const std::vector<Real>& x = anti.abscissa();
   const std::vector<Real>& w = anti.weights();
   CHECK_EQUAL(x.size(), static_cast<std::size_t>(N + 1));
   CHECK_LE(Real(-1), x[0]);
   CHECK_LE(x[N], Real(1));
   Real sum = 0;
   for (std::size_t i = 0; i < x.size(); ++i)
   {
      CHECK_LE(Real(0), w[i]);
      sum += w[i];
   }
   CHECK_ULP_CLOSE(pow(Real(2), alpha + beta + 1) * boost::math::beta(alpha + 1, beta + 1), sum, 16);
   for (std::size_t i = 0; i < N; ++i)
   {
      CHECK_LE(x[i], rule.abscissa()[i]);
      CHECK_LE(rule.abscissa()[i], x[i + 1]);
   }
   // p(x) = sum c_k (1 + x)^k with arbitrary c_k, k <= 2N + 1, whose integral is
   // sum c_k 2^(alpha + beta + k + 1) B(alpha + 1, beta + k + 1).
   std::vector<Real> c(2 * N + 2);
   Real exact = 0;
   for (unsigned k = 0; k < c.size(); ++k)
   {
      c[k] = Real((k * 7919) % 17) / 8 - 1;
      exact += c[k] * pow(Real(2), alpha + beta + k + 1) * boost::math::beta(alpha + 1, beta + k + 1);
   }
   auto p = [&c](const Real& t) { Real r = 0; for (std::size_t k = c.size(); k-- > 0;) r = r * (1 + t) + c[k]; return r; };
   Real L1;
   const Real gauss = rule.integrate(p, nullptr, &L1);
   Real anti_result = 0;
   for (std::size_t i = 0; i < x.size(); ++i)
   {
      anti_result += w[i] * p(x[i]);
      L1 += w[i] * abs(p(x[i]));
   }
   CHECK_ABSOLUTE_ERROR(Real(0), (exact - gauss) + (exact - anti_result), (N + 8) * std::numeric_limits<Real>::epsilon() * L1);
}

// Laurie 1996, Theorem 3: the anti-Gauss rule is internal iff both of these are nonnegative. Zero puts a node on an
// endpoint, which the rule rejects too.
template <class Real>
void test_anti_gauss_internality()
{
   using family = boost::math::detail::jacobi_family<Real>;
   using polynomial = boost::math::detail::orthogonal_polynomial<Real, family>;
   const Real parameters[] = { -0.75, -0.55, -0.5, -0.25, 0, 0.5, 3 };
   for (Real a : parameters)
   {
      for (Real b : parameters)
      {
         for (unsigned n = 1; n <= 40; ++n)
         {
            const Real right = (2 * a + 1) * n * n + (2 * a + 1) * (a + b + 1) * n + (a + 1) * (a + b) * (a + b + 1) / 2;
            const Real left = (2 * b + 1) * n * n + (2 * b + 1) * (a + b + 1) * n + (b + 1) * (a + b) * (a + b + 1) / 2;
            const family f(n + 1, a, b);
            CHECK_EQUAL(polynomial::anti_gauss_is_internal(n + 1, f, Real(1)), right > 0);
            CHECK_EQUAL(polynomial::anti_gauss_is_internal(n + 1, f, Real(-1)), left > 0);
         }
      }
   }
   // Chebyshev weights of the first kind put anti-Gauss nodes exactly on +-1.
   auto internal = [](unsigned n, Real a, Real b)
   {
      const family f(n + 1, a, b);
      return polynomial::anti_gauss_is_internal(n + 1, f, Real(1)) && polynomial::anti_gauss_is_internal(n + 1, f, Real(-1));
   };
   for (unsigned n = 1; n <= 200; n += 7)
   {
      CHECK_EQUAL(internal(n, Real(-0.5), Real(-0.5)), false);
      CHECK_EQUAL(internal(n, Real(0), Real(0)), true);
      CHECK_EQUAL(internal(n, Real(0.5), Real(0.5)), true);
      CHECK_EQUAL(internal(n, Real(-0.75), Real(0)), false);
   }
   // Without an internal anti-Gauss rule there is no estimate, but still an integral.
   using ignore_policy = boost::math::policies::policy<boost::math::policies::evaluation_error<boost::math::policies::ignore_error> >;
   gauss_jacobi<Real, ignore_policy> chebyshev(10, Real(-0.5), Real(-0.5));
   auto f = [](const Real& x) { return x * x; };
   Real error = 0;
   CHECK_ULP_CLOSE(boost::math::constants::half_pi<Real>(), chebyshev.integrate(f, &error, nullptr), 8);
   CHECK_NAN(error);
#ifndef BOOST_NO_EXCEPTIONS
   gauss_jacobi<Real> throwing(10, Real(-0.5), Real(-0.5));
   CHECK_THROW(throwing.integrate(f, &error, nullptr), boost::math::evaluation_error);
   CHECK_ULP_CLOSE(boost::math::constants::half_pi<Real>(), throwing.integrate(f, nullptr, nullptr), 8);
#endif
}

// For smooth f the estimate is close to the true error until that reaches rounding.
template <class Real>
void test_error_estimate(unsigned N)
{
   using std::abs;
   using std::cos;
   using std::exp;
   struct example
   {
      Real alpha;
      Real beta;
      Real (*f)(const Real&);
      Real expected;
   };
   const example examples[] = {
      { Real(0.5), Real(0.5), [](const Real& x) { return cos(x); }, static_cast<Real>(1.382459687384168525766283327693941037437L) },
      { Real(1), Real(1), [](const Real& x) { return cos(x); }, static_cast<Real>(1.204674715759027157006262856749289583561L) },
      { Real(0.3), Real(-0.4), [](const Real& x) { return exp(x); }, static_cast<Real>(2.113791650902447106944400891924615994062L) },
      { Real(0), Real(0), [](const Real& x) { return 1 / (2 - x); }, static_cast<Real>(1.098612288668109691395245236922525704647L) },
   };
   for (const example& e : examples)
   {
      gauss_jacobi<Real> rule(N, e.alpha, e.beta);
      Real error, L1, L1_alone;
      const Real Q = rule.integrate(e.f, &error, &L1);
      CHECK_ULP_CLOSE(Q, rule.integrate(e.f, nullptr, &L1_alone), 2);
      CHECK_ULP_CLOSE(L1, L1_alone, 2);
      const Real true_error = abs(Q - e.expected);
      if (true_error > 1000 * std::numeric_limits<Real>::epsilon() * L1)
      {
         CHECK_LE(error, 2 * true_error);
         CHECK_LE(true_error, 2 * error);
      }
      CHECK_LE(2 * std::numeric_limits<Real>::epsilon() * L1, error);
   }
}

// No higher precision because it takes a long time because it has to
// compute the nodes/weights on the fly!
int main()
{
   test_moments<float, 6>(8 * std::numeric_limits<float>::epsilon(), "float");
   test_moments<double, 6>(5 * std::numeric_limits<double>::epsilon(), "double");
   test_moments<double, 10>(5 * std::numeric_limits<double>::epsilon(), "double");

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_moments<long double, 6>(4 * std::numeric_limits<double>::epsilon(), "long double");
#endif
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_moments<boost::math::concepts::real_concept, 6>(4e-16, "real concept");
#endif


   test_complex<std::complex<double>, 30>(4 * std::numeric_limits<double>::epsilon());

   test_vector_valued<30>();
   test_vector_valued<40>();

   for (unsigned n : { 1u, 2u, 5u, 10u, 15u })
   {
      test_anti_gauss<double>(n, 0, 0);
      test_anti_gauss<double>(n, 0.5, 0.5);
      test_anti_gauss<double>(n, 0.3, -0.4);
      test_anti_gauss<double>(n, 2.5, 1);
   }
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_anti_gauss<long double>(9, 0.25, -0.25);
#endif
   test_anti_gauss_internality<double>();
   for (unsigned n : { 2u, 3u, 5u, 8u, 12u, 20u })
   {
      test_error_estimate<double>(n);
   }
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_error_estimate<long double>(6);
#endif

   return boost::math::test::report_errors();
}
