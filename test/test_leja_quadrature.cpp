// Copyright Nick Thompson, 2026
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#include "math_unit_test.hpp"
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <set>
#include <stdexcept>
#include <thread>
#include <vector>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/leja.hpp>
#include <boost/math/quadrature/sparse_grid.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/complex128.hpp>
#endif

using boost::math::quadrature::leja_hermite;
using boost::math::quadrature::leja_jacobi;
using boost::math::quadrature::leja_laguerre;
using boost::math::quadrature::leja_legendre;

// The n-point rule, sum_i w_i f(x_i), and L1 = sum_i |w_i| |f(x_i)|, for a real or complex f.
template <class Rule, class F>
auto n_point(const Rule& rule, F f, unsigned n, typename Rule::value_type* L1 = nullptr)
{
   using std::abs;
   const std::vector<typename Rule::value_type>& w = rule.weights(n);
   decltype(f(w[0])) Q(0);
   typename Rule::value_type l = 0;
   for (unsigned i = 0; i < n; ++i)
   {
      const auto y = f(rule.abscissa()[i]);
      Q += w[i] * y;
      l += abs(w[i]) * abs(y);
   }
   if (L1)
      *L1 = l;
   return Q;
}

// The integral of x^k over [-1, 1], or of x^k exp(-x^2) over the real line, Gamma((k + 1)/2) for even k.
template <class Real>
Real moment(const leja_legendre<Real>&, unsigned k)
{
   return k % 2 ? Real(0) : Real(2) / (k + 1);
}

template <class Real>
Real moment(const leja_hermite<Real>&, unsigned k)
{
   Real m = boost::math::constants::root_pi<Real>();
   for (unsigned j = 2; j <= k; j += 2)
      m *= Real(j - 1) / 2;
   return k % 2 ? Real(0) : m;
}


// Rounding relative to L1, in units of eps: the small Laguerre weights far out, where x^k is large, are the least accurate.
template <class Rule>
unsigned slack(const Rule&, unsigned n)
{
   return 8 * n;
}

template <class Real>
unsigned slack(const leja_laguerre<Real>&, unsigned n)
{
   return 32 * n;
}

// The n-point rule is exact for polynomials of degree n - 1, up to rounding relative to L1 = sum |w_i x_i^k|; moment(k) is
// the integral of x^k. The sum of |w_i| is at most mass times the integral of 1.
template <class Rule, class Moment>
void test_exactness(const Rule& rule, Moment moment, unsigned N, double mass = 1.25)
{
   using Real = typename Rule::value_type;
   using std::abs;
   using std::pow;
   const Real eps = std::numeric_limits<Real>::epsilon();
   for (unsigned n = 1; n <= N; ++n)
   {
      CHECK_EQUAL(rule.weights(n).size(), std::size_t(n));
      Real sum = 0;
      for (const Real& w : rule.weights(n))
         sum += abs(w);
      // The weights can be negative, but not by much.
      CHECK_LE(sum, Real(mass) * moment(0));
      for (unsigned k = 0; k < n; ++k)
      {
         Real L1;
         const Real Q = n_point(rule, [k](const Real& x) { return pow(x, Real(k)); }, n, &L1);
         CHECK_LE(abs(Q - moment(k)), slack(rule, n) * eps * L1);
      }
   }
}

// A rule computed to 40 points at once is the rule computed one point at a time, and its copy.
template <class Rule>
void check_nested(const Rule& rule)
{
   const Rule steps(rule);
   rule.weights(40);
   const Rule copy(rule);
   CHECK_EQUAL(copy.abscissa().size(), std::size_t(40));
   // Assigned, a rule keeps the points and extends them as the original does.
   Rule assigned(steps);
   assigned = copy;
   CHECK_EQUAL(assigned.weights(41)[40], rule.weights(41)[40]);
   for (unsigned n = 1; n <= 40; ++n)
   {
      const std::vector<typename Rule::value_type>& w = steps.weights(n);
      CHECK_EQUAL(steps.abscissa()[n - 1], rule.abscissa()[n - 1]);
      for (unsigned i = 0; i < n; ++i)
      {
         CHECK_EQUAL(w[i], rule.weights(n)[i]);
         CHECK_EQUAL(copy.weights(n)[i], rule.weights(n)[i]);
      }
   }
}

template <class Real>
void test_exactness(unsigned n = std::is_same<Real, float>::value ? 16 : 32)
{
   using boost::math::tgamma;
   const leja_legendre<Real> legendre;
   const leja_hermite<Real> hermite;
   test_exactness(legendre, [&](unsigned k) { return moment(legendre, k); }, n);
   test_exactness(hermite, [&](unsigned k) { return moment(hermite, k); }, n);
   for (Real a : {Real(0), Real(-0.5L), Real(2.5L)})
      test_exactness(leja_laguerre<Real>(a), [a](unsigned k) { return tgamma(k + a + 1); }, n);
   // n-point Gauss-Jacobi integrates x^k, k < 2n, to rounding.
   for (Real a : {Real(-0.5L), Real(0.5L)})
   {
      for (Real b : {Real(-0.75L), Real(2)})
      {
         boost::math::quadrature::gauss_jacobi<Real> gauss(n, a, b);
         // A negative exponent places the points as if it were 0, and the weights are more negative.
         test_exactness(leja_jacobi<Real>(a, b), [&](unsigned k) { return gauss.integrate([k](const Real& x) { using std::pow; return pow(x, Real(k)); }); }, n, 1.5);
      }
   }
}

template <class Real>
void test_nested()
{
   leja_legendre<Real> large;
   leja_hermite<Real> h_large;
   check_nested(large);
   check_nested(h_large);
   check_nested(leja_laguerre<Real>(Real(1.5)));
   check_nested(leja_jacobi<Real>(Real(-0.5), Real(2)));
   // The first points: the mean of the weight, then 0, -1, 1 for Legendre, where exact ties go to the left.
   CHECK_EQUAL(large.abscissa()[0], Real(0));
   CHECK_ULP_CLOSE(-1 / boost::math::constants::root_three<Real>(), large.abscissa()[3], 2);
   CHECK_EQUAL(h_large.abscissa()[0], Real(0));
   CHECK_ULP_CLOSE(Real(-1), h_large.abscissa()[1], 2);
   leja_laguerre<Real> laguerre;
   laguerre.weights(2);
   CHECK_EQUAL(laguerre.abscissa()[1], Real(0));
   leja_jacobi<Real> chebyshev(Real(-0.5), Real(-0.5));
   chebyshev.weights(3);
   CHECK_EQUAL(chebyshev.abscissa()[1], Real(-1));
   CHECK_EQUAL(chebyshev.abscissa()[2], Real(1));
}

// Every precision picks the same point at every step.
template <class Real, template <class, class> class Rule, class... Args>
void test_same_points(unsigned points, const Args&... args)
{
   using std::abs;
   Rule<Real, boost::math::policies::policy<>> rule{Real(args)...};
   Rule<double, boost::math::policies::policy<>> reference{double(args)...};
   rule.weights(points);
   reference.weights(points);
   const Real tol = 8 * (std::max)(Real(std::numeric_limits<double>::epsilon()), std::numeric_limits<Real>::epsilon());
   for (unsigned i = 0; i < points; ++i)
      CHECK_ABSOLUTE_ERROR(Real(reference.abscissa()[i]), rule.abscissa()[i], tol * (std::max)(Real(1), abs(rule.abscissa()[i])));
}

template <class Real>
void test_smooth(unsigned n)
{
   using std::cos;
   using std::exp;
   using std::sin;
   const Real eps = std::numeric_limits<Real>::epsilon();
   auto f = [](const Real& x) { return cos(3 * x); };
   leja_legendre<Real> rule;
   CHECK_ABSOLUTE_ERROR(2 * sin(Real(3)) / 3, n_point(rule, f, n), 16 * eps);
   // Hermite needs twice as many points: Leja rules are exact to degree n - 1, not 2n - 1.
   leja_hermite<Real> h_rule;
   CHECK_ABSOLUTE_ERROR(boost::math::constants::root_pi<Real>() * exp(Real(-9) / 4), n_point(h_rule, f, 2 * n), 16 * eps);
}

template <class Real>
void test_adaptive()
{
   using std::abs;
   using std::atan;
   using std::cos;
   using std::exp;
   using std::sin;
   using std::sqrt;
   const Real eps = std::numeric_limits<Real>::epsilon();
   leja_legendre<Real> rule;
   leja_hermite<Real> h_rule;
   leja_laguerre<Real> l_rule;
   leja_laguerre<Real> l_half(Real(-0.5L));
   leja_jacobi<Real> chebyshev(Real(-0.5L), Real(-0.5L));
   Real error, L1;
   // On these integrands the error estimate bounds the error, or the error is at the level of rounding.
   auto check = [&](const Real& expected, const Real& Q, const Real& tolerance) {
      CHECK_LE(error, tolerance * L1);
      CHECK_LE(abs(Q - expected), (std::max)(error, 64 * eps * L1));
   };
   Real Q = rule.integrate([](const Real& x) { return cos(3 * x); }, Real(1e-12L), &error, &L1);
   check(2 * sin(Real(3)) / 3, Q, Real(1e-12L));
   Q = rule.integrate([](const Real& x) { return 1 / (1 + 16 * x * x); }, Real(1e-6L), &error, &L1);
   check(atan(Real(4)) / 2, Q, Real(1e-6L));
   Q = rule.integrate([](const Real& x) { return exp(x); }, boost::math::tools::root_epsilon<Real>(), &error, &L1);
   check(exp(Real(1)) - exp(Real(-1)), Q, boost::math::tools::root_epsilon<Real>());
   Q = h_rule.integrate([](const Real& x) { return cos(3 * x); }, Real(1e-12L), &error, &L1);
   check(boost::math::constants::root_pi<Real>() * exp(Real(-9) / 4), Q, Real(1e-12L));
   Q = h_rule.integrate([](const Real& x) { return 1 / (1 + x * x / 4); }, Real(1e-8L), &error, &L1);
   // 2 pi e^4 erfc(2)
   check(Real(1.604698360911363269322784609663024807L), Q, Real(1e-8L));
   Q = l_rule.integrate([](const Real& x) { return cos(x); }, Real(1e-10L), &error, &L1);
   check(Real(0.5L), Q, Real(1e-10L));
   Q = l_rule.integrate([](const Real& x) { return 1 / (1 + x); }, Real(1e-6L), &error, &L1);
   // e E_1(1)
   check(Real(0.5963473623231940743410784993692793760742L), Q, Real(1e-6L));
   // Gamma(1/2) cos(pi/8) / 2^(1/4) and pi J_0(1).
   Q = l_half.integrate([](const Real& x) { return cos(x); }, Real(1e-10L), &error, &L1);
   check(Real(1.376996331853153438664376624047155467474L), Q, Real(1e-10L));
   Q = chebyshev.integrate([](const Real& x) { return cos(x); }, Real(1e-12L), &error, &L1);
   check(Real(2.403939430634412998273324892591511223785L), Q, Real(1e-12L));
}

// e^{ix} over [-1, 1] with n points, and with exp(-x^2) with h_n.
template <class Real, class Complex = std::complex<Real>>
void test_complex(unsigned n, unsigned h_n)
{
   using std::cos;
   using std::exp;
   using std::sin;
   const Real eps = std::numeric_limits<Real>::epsilon();
   const Real tolerance = boost::math::tools::root_epsilon<Real>();
   auto f = [](const Real& x) { return Complex(cos(x), sin(x)); };
   leja_legendre<Real> rule;
   leja_hermite<Real> h_rule;
   Real error, L1;
   Complex Q = n_point(rule, f, n);
   CHECK_ABSOLUTE_ERROR(2 * sin(Real(1)), Real(Q.real()), 8 * eps);
   CHECK_ABSOLUTE_ERROR(Real(0), Real(Q.imag()), 8 * eps);
   Q = rule.integrate(f, tolerance, &error, &L1);
   CHECK_ABSOLUTE_ERROR(2 * sin(Real(1)), Real(Q.real()), (std::max)(error, 8 * eps));
   CHECK_LE(error, tolerance * L1);
   Q = n_point(h_rule, f, h_n);
   const Real expected = boost::math::constants::root_pi<Real>() * exp(Real(-1) / 4);
   CHECK_ABSOLUTE_ERROR(expected, Real(Q.real()), 8 * eps);
   CHECK_ABSOLUTE_ERROR(Real(0), Real(Q.imag()), 8 * eps);
   Q = h_rule.integrate(f, tolerance, &error, &L1);
   CHECK_ABSOLUTE_ERROR(expected, Real(Q.real()), (std::max)(error, 8 * eps));
   CHECK_LE(error, tolerance * L1);
}

// A vector-valued integrand is integrated componentwise, and the adaptive rule waits for the slowest component.
template <class Real, class Complex = std::complex<Real>>
void test_vector()
{
   using std::abs;
   using std::cos;
   using std::exp;
   using std::sin;
   const Real eps = std::numeric_limits<Real>::epsilon();
   const Real tolerance = Real(1e-6L);
   leja_legendre<Real> rule;
   unsigned calls = 0;
   auto f = [&calls](const Real& x) {
      ++calls;
      return std::array<Real, 3>{x * x, cos(3 * x), 1 / (1 + x * x)};
   };
   const std::array<Real, 3> expected{Real(2) / 3, 2 * sin(Real(3)) / 3, boost::math::constants::half_pi<Real>()};
   Real error, L1;
   // A constant converges at 4 points, and L1 is the max norm 2 times the sum of |w_i|.
   std::array<Real, 2> C = rule.integrate([](const Real&) { return std::array<Real, 2>{Real(1), Real(-2)}; }, tolerance, &error, &L1);
   CHECK_ULP_CLOSE(Real(-4), C[1], 4);
   CHECK_ULP_CLOSE(2 * n_point(rule, [](const Real&) { return Real(1); }, 4), L1, 4);
   std::array<Real, 3> Q = rule.integrate(f, tolerance, &error, &L1);
   CHECK_LE(error, tolerance * L1);
   for (std::size_t j = 0; j < 3; ++j)
      CHECK_LE(abs(Q[j] - expected[j]), (std::max)(error, 64 * eps * L1));
   // A NaN component is not hidden by the max norm.
   rule.integrate([](const Real&) { return std::array<Real, 2>{std::numeric_limits<Real>::quiet_NaN(), Real(1)}; }, tolerance, &error, &L1);
   CHECK_NAN(L1);
   unsigned slowest = 0;
   rule.integrate([&slowest](const Real& x) { ++slowest; return 1 / (1 + x * x); }, tolerance);
   CHECK_LE(slowest, calls);

   // A std::vector of complex values, with the Hermite weight.
   leja_hermite<Real> h_rule;
   auto g = [](const Real& x) { return std::vector<Complex>{Complex(cos(x), sin(x)), Complex(x * x, Real(0))}; };
   std::vector<Complex> P = h_rule.integrate(g, tolerance, &error, &L1);
   CHECK_EQUAL(P.size(), std::size_t(2));
   CHECK_LE(error, tolerance * L1);
   CHECK_ABSOLUTE_ERROR(boost::math::constants::root_pi<Real>() * exp(Real(-1) / 4), Real(P[0].real()), (std::max)(error, 8 * eps));
   CHECK_ABSOLUTE_ERROR(Real(0), Real(P[0].imag()), (std::max)(error, 8 * eps));
   CHECK_ABSOLUTE_ERROR(boost::math::constants::root_pi<Real>() / 2, Real(P[1].real()), 8 * eps);
}

// The sparse grid, on integrands with known integrals. It stops when its error estimate is below tolerance * L1.
template <class Real, class Complex = std::complex<Real>>
void test_sparse_grid(const Real& tolerance)
{
   using boost::math::quadrature::sparse_grid;
   using std::abs;
   using std::cos;
   using std::exp;
   using std::pow;
   using std::sin;
   using point = std::vector<Real>;
   const Real eps = std::numeric_limits<Real>::epsilon();
   const Real pi = boost::math::constants::pi<Real>();
   leja_legendre<Real> rule;
   leja_hermite<Real> h_rule;
   Real error, L1;
   auto check = [&](const Real& expected, const Real& Q) {
      CHECK_LE(error, tolerance * L1);
      CHECK_LE(abs(Q - expected), (std::max)(error, 64 * eps * L1));
   };

   // Constants, odd functions, and x^2, whose rules with 1 and 2 points agree.
   check(Real(32), sparse_grid([](const point&) { return Real(1); }, rule, 5, 1000, tolerance, &error, &L1));
   check(pow(pi, Real(5) / 2), sparse_grid([](const point&) { return Real(1); }, h_rule, 5, 1000, tolerance, &error, &L1));
   check(Real(0), sparse_grid([](const point& x) { return x[0] * x[1] * x[1] * x[1]; }, rule, 3, 1000, tolerance, &error, &L1));
   check(Real(2) / 3, sparse_grid([](const point& x) { return x[0] * x[0]; }, rule, 1, 1000, tolerance, &error, &L1));
   check(pi / 2, sparse_grid([](const point& x) { return x[1] * x[1]; }, h_rule, 2, 1000, tolerance, &error, &L1));
   // With the Laguerre weight, exp(-x_0/2 - x_1/3 - x_2/4) integrates to 2/3 * 3/4 * 4/5, and x_0 x_1, which vanishes on
   // the axes through 0, to 1.
   check(Real(2) / 5, sparse_grid([](const point& x) { return exp(-x[0] / 2 - x[1] / 3 - x[2] / 4); }, leja_laguerre<Real>(), 3,
                                  100000, tolerance, &error, &L1));
   check(Real(1), sparse_grid([](const point& x) { return x[0] * x[1]; }, leja_laguerre<Real>(), 2, 1000, tolerance, &error, &L1));
   // f may take a std::vector, be generic, or take a std::array, as a lambda or a function pointer.
   static const auto e = [](const auto& x) { return exp(x[0] / 2 - x[1] + x[2] / 4); };
   const Real Qv = sparse_grid([](const point& x) { return e(x); }, rule, 3, 10000, tolerance);
   CHECK_EQUAL(Qv, sparse_grid(e, rule, 3, 10000, tolerance));
   CHECK_EQUAL(Qv, sparse_grid([](const std::array<Real, 3>& x) { return e(x); }, rule, 3, 10000, tolerance));
   CHECK_EQUAL(Qv, sparse_grid(+[](std::array<Real, 3> x) { return e(x); }, rule, 3, 10000, tolerance));
   // With the Chebyshev weight, cos(x_0 + x_1/2) integrates to pi^2 J_0(1) J_0(1/2), the same from a rule that grows
   // as from one that has its points already.
   auto c2 = [](const point& x) { return cos(x[0] + x[1] / 2); };
   const leja_jacobi<Real> chebyshev(Real(-0.5L), Real(-0.5L));
   const Real Qc = sparse_grid(c2, chebyshev, 2, 10000, tolerance, &error, &L1);
   check(Real(7.087510228266830598984765096456477556883L), Qc);
   chebyshev.weights(63);
   CHECK_EQUAL(Qc, sparse_grid(c2, chebyshev, 2, 10000, tolerance, &error, &L1));
   // A NaN stops the grid after the refinement that found it.
   unsigned nan_calls = 0;
   CHECK_NAN(sparse_grid([&](const point& x) { ++nan_calls; return x[0] > 0 ? Real(1) : std::numeric_limits<Real>::quiet_NaN(); },
                         rule, 3, 1000, tolerance, &error, &L1));
   CHECK_LE(nan_calls, 1u);
   // Out of calls before the check is done: a NaN error estimate.
   sparse_grid([](const point& x) { return x[0] * x[0]; }, rule, 1, 1, tolerance, &error, &L1);
   CHECK_NAN(error);
   // A product of polynomials of degree 4, exact once every level is 3.
   auto p = [](const Real& t) { return 1 + t - 2 * t * t + t * t * t * t; };
   check(Real(16 * 16 * 16) / (15 * 15 * 15),
         sparse_grid([&](const point& x) { return p(x[0]) * p(x[1]) * p(x[2]); }, rule, 3, 1000, tolerance, &error, &L1));

   // Smooth and anisotropic: exp(sum_i x_i / 2^i) and, with the Hermite weight, cos(sum_i x_i / (i + 1)).
   const unsigned d = 6;
   Real expected = 1;
   for (unsigned i = 0; i < d; ++i)
      expected *= (exp(Real(1) / (1 << i)) - exp(-Real(1) / (1 << i))) * (1 << i);
   std::set<point> seen;
   unsigned calls = 0;
   auto f = [&](const point& x) {
      ++calls;
      CHECK_EQUAL(seen.insert(x).second, true);
      Real s = 0;
      for (unsigned i = 0; i < d; ++i)
         s += x[i] / (1 << i);
      return exp(s);
   };
   check(expected, sparse_grid(f, rule, d, 100000, tolerance, &error, &L1));
   Real a2 = 0;
   for (unsigned i = 0; i < 4; ++i)
      a2 += 1 / Real((i + 1) * (i + 1));
   auto g = [](const point& x) {
      Real s = 0;
      for (unsigned i = 0; i < 4; ++i)
         s += x[i] / (i + 1);
      return Complex(cos(s), sin(s));
   };
   Complex Q = sparse_grid(g, h_rule, 4, 100000, tolerance, &error, &L1);
   check(pi * pi * exp(-a2 / 4), Real(Q.real()));
   CHECK_LE(Real(abs(Q.imag())), (std::max)(error, 64 * eps * L1));
   // A vector is integrated componentwise, with the max norm.
   auto h = [&](const point& x) { return std::array<Complex, 2>{g(x), Complex(x[0] * x[0], Real(0))}; };
   std::array<Complex, 2> R = sparse_grid(h, h_rule, 4, 100000, tolerance, &error, &L1);
   check(pi * pi * exp(-a2 / 4), Real(R[0].real()));
   check(pi * pi / 2, Real(R[1].real()));

   // Out of calls: never more than max_calls, and an error estimate above the tolerance.
   calls = 0;
   seen.clear();
   sparse_grid(f, rule, d, 50, tolerance, &error, &L1);
   CHECK_LE(calls, 50u);
   CHECK_LE(tolerance * L1, error);
}

// Functions that vanish at the first point, the mean of the weight, and on the axes through it.
void test_vanishing()
{
   using boost::math::quadrature::sparse_grid;
   using point = std::vector<double>;
   const double eps = std::numeric_limits<double>::epsilon();
   const double tolerance = boost::math::tools::root_epsilon<double>();
   const double e = std::exp(1.0);
   double error, L1;
   auto check = [&](const double& expected, const double& Q) {
      CHECK_LE(error, tolerance * L1);
      CHECK_LE(std::abs(Q - expected), (std::max)(error, 64 * eps * L1));
   };
   leja_legendre<double> rule;
   check(4.0 / 9, sparse_grid([](const point& x) { return x[0] * x[0] * x[1] * x[1]; }, rule, 2, 100000, tolerance, &error, &L1));
   check((2 * std::sin(1.0) - 2) * (2 * std::sin(1.0) - 2),
         sparse_grid([](const point& x) { return (std::cos(x[0]) - 1) * (std::cos(x[1]) - 1); }, rule, 2, 100000, tolerance, &error, &L1));
   check(4 / (e * e), sparse_grid([](const point& x) { return x[0] * x[1] * std::exp(x[0] + x[1]); }, rule, 2, 100000, tolerance, &error, &L1));
   check(8 / (e * e * e), sparse_grid([](const point& x) { return x[0] * x[1] * x[2] * std::exp(x[0] + x[1] + x[2]); }, rule, 3, 100000,
                                      tolerance, &error, &L1));
   check(boost::math::constants::pi<double>() / 4,
         sparse_grid([](const point& x) { return x[0] * x[0] * x[1] * x[1]; }, leja_hermite<double>(), 2, 100000, tolerance, &error, &L1));
   // 0 everywhere: 0 after dimension + 1 rounds, which reach the product x_0^2 ... x_3^2.
   for (unsigned d = 1; d <= 4; ++d)
   {
      unsigned calls = 0;
      CHECK_EQUAL(sparse_grid([&](const point&) { ++calls; return 0.0; }, rule, d, 100000, tolerance, &error, &L1), 0.0);
      CHECK_LE(calls, 1000u);
      CHECK_EQUAL(error, 0.0);
      CHECK_EQUAL(L1, 0.0);
   }
   check(16.0 / 81, sparse_grid([](const point& x) { return x[0] * x[0] * x[1] * x[1] * x[2] * x[2] * x[3] * x[3]; }, rule, 4, 100000,
                                tolerance, &error, &L1));
}

// Whether integrate(f, a, b) is declared.
template <class Rule, class = void>
struct has_bounds : std::false_type {};
template <class Rule>
struct has_bounds<Rule, std::void_t<decltype(std::declval<const Rule&>().integrate(std::declval<double (*)(double)>(), 0.0, 2.0))>> : std::true_type {};
static_assert(has_bounds<leja_legendre<double>>::value && has_bounds<leja_jacobi<double>>::value, "");
static_assert(!has_bounds<leja_hermite<double>>::value && !has_bounds<leja_laguerre<double>>::value, "");

// integrate(f, a, b) and sparse_grid(f, rule, a, b) map [-1, 1] affinely onto [a, b].
template <class Real, class Complex = std::complex<Real>>
void test_bounds(const Real& tolerance)
{
   using boost::math::quadrature::sparse_grid;
   using std::abs;
   using std::cos;
   using std::exp;
   using std::pow;
   using std::sin;
   const Real eps = std::numeric_limits<Real>::epsilon();
   const leja_legendre<Real> rule;
   Real error, L1, error2, L2;
   auto check = [&](const Real& expected, const Real& Q) {
      CHECK_LE(error, tolerance * L1);
      CHECK_LE(abs(Q - expected), (std::max)(error, 64 * eps * L1));
   };
   // Every form resolves to the intended overload.
   auto e = [](const Real& x) { return exp(x); };
   const Real e1 = exp(Real(1)) - exp(Real(-1)), e02 = exp(Real(2)) - 1;
   const Real bound = (std::max)(boost::math::tools::root_epsilon<Real>(), Real(1e-8));
   CHECK_LE(abs(rule.integrate(e) - e1), bound);
   CHECK_LE(abs(rule.integrate(e, 1e-10) - e1), bound);
   CHECK_LE(abs(rule.integrate(e, 1e-10, &error, &L1) - e1), bound);
   CHECK_LE(abs(rule.integrate(e, 0.0, 2.0) - e02), bound * e02);
   CHECK_LE(abs(rule.integrate(e, 0.0, 2.0, 1e-10) - e02), bound * e02);
   CHECK_LE(abs(rule.integrate(e, 0.0, 2.0, 1e-10, &error, &L1) - e02), bound * e02);
   CHECK_LE(abs(rule.integrate(e, 0, 2) - e02), bound * e02);

   // A polynomial of degree 4 on [0, 3], exact from 5 points, and exp on [-2, 5].
   auto p = [](const Real& x) { return 1 + x - 2 * x * x + x * x * x * x; };
   check(Real(381) / 10, rule.integrate(p, Real(0), Real(3), tolerance, &error, &L1));
   const Real Q = rule.integrate(e, Real(-2), Real(5), tolerance, &error, &L1);
   check(exp(Real(5)) - exp(Real(-2)), Q);
   // The error and L1 of the rule for f(x(t)) on [-1, 1], times |b - a|/2.
   const Real Qt = rule.integrate([](const Real& t) { return exp(Real(1.5) + Real(3.5) * t); }, tolerance, &error2, &L2);
   CHECK_ULP_CLOSE(Real(3.5) * Qt, Q, 4);
   CHECK_ULP_CLOSE(Real(3.5) * error2, error, 4);
   CHECK_ULP_CLOSE(Real(3.5) * L2, L1, 4);
   // a > b: minus the integral over [b, a], with a positive error and L1.
   const Real R = rule.integrate(e, Real(5), Real(-2), tolerance, &error2, &L2);
   CHECK_LE(error2, tolerance * L2);
   CHECK_LE(abs(R + exp(Real(5)) - exp(Real(-2))), (std::max)(error2, 64 * eps * L2));
   // a == b: 0, of the shape of f(a), from one call.
   unsigned calls = 0;
   CHECK_EQUAL(rule.integrate([&](const Real& x) { ++calls; return x; }, Real(1.5), Real(1.5), tolerance, &error, &L1), Real(0));
   CHECK_EQUAL(calls, 1u);
   CHECK_EQUAL(error, Real(0));
   CHECK_EQUAL(L1, Real(0));
   const std::vector<Real> Z = rule.integrate([](const Real& x) { return std::vector<Real>{x, x * x}; }, Real(1), Real(1));
   CHECK_EQUAL(Z.size(), std::size_t(2));
   CHECK_EQUAL(Z[0], Real(0));
   CHECK_EQUAL(Z[1], Real(0));
   const Real inf = std::numeric_limits<Real>::infinity();
   CHECK_THROW(rule.integrate(e, Real(0), inf), std::domain_error);
   CHECK_THROW(rule.integrate(e, -inf, Real(0)), std::domain_error);
   CHECK_THROW(rule.integrate(e, std::numeric_limits<Real>::quiet_NaN(), Real(0)), std::domain_error);

   // Complex and vector f: e^{ix} over [0, pi] is 2i, and (x, x^2) over [0, 2] is (2, 8/3).
   const Real pi = boost::math::constants::pi<Real>();
   const Complex C = rule.integrate([](const Real& x) { return Complex(cos(x), sin(x)); }, Real(0), pi, tolerance, &error, &L1);
   CHECK_LE(Real(abs(C - Complex(0, 2))), (std::max)(error, 64 * eps * L1));
   const std::array<Real, 2> V = rule.integrate([](const Real& x) { return std::array<Real, 2>{x, x * x}; }, Real(0), Real(2), tolerance, &error, &L1);
   check(Real(2), V[0]);
   check(Real(8) / 3, V[1]);

   // The integral of x^k (4 - x)^alpha (x - 1)^beta over [1, 4] is 3^(alpha + beta + 1) B(alpha + 1, beta + 1) for k = 0,
   // and that times 1 + 3 (beta + 1)/(alpha + beta + 2) for k = 1. Over [4, 1], alpha goes with 1, and the sign flips.
   for (Real alpha : {Real(0.5), Real(2)})
   {
      for (Real beta : {Real(-0.25), Real(1.5)})
      {
         const leja_jacobi<Real> jacobi(alpha, beta);
         const Real m = pow(Real(3), alpha + beta + 1) * boost::math::beta(alpha + 1, beta + 1);
         const Real mean = 1 + 3 * (beta + 1) / (alpha + beta + 2);
         check(m, jacobi.integrate([](const Real&) { return Real(1); }, Real(1), Real(4), tolerance, &error, &L1));
         check(m * mean, jacobi.integrate([](const Real& x) { return x; }, Real(1), Real(4), tolerance, &error, &L1));
         check(-m * (5 - mean), jacobi.integrate([](const Real& x) { return x; }, Real(4), Real(1), tolerance, &error, &L1));
         // The error and L1 scale by (3/2)^(alpha + beta + 1).
         const Real Qj = jacobi.integrate(e, Real(1), Real(4), tolerance, &error, &L1);
         const Real s = pow(Real(1.5), alpha + beta + 1);
         CHECK_ULP_CLOSE(s * jacobi.integrate([](const Real& t) { return exp(Real(2.5) + Real(1.5) * t); }, tolerance, &error2, &L2), Qj, 8);
         CHECK_ULP_CLOSE(s * error2, error, 8);
         CHECK_ULP_CLOSE(s * L2, L1, 8);
         // x_0 in a sparse grid on [1, 4]^2.
         check(m * m * mean, sparse_grid([](const std::vector<Real>& x) { return x[0]; }, jacobi, std::vector<Real>{1, 1},
                                         std::vector<Real>{4, 4}, 10000, tolerance, &error, &L1));
      }
   }

   // exp(x_0/2 - x_1 + x_2/4) on [0, 1] x [-1, 2] x [3, 5] with std::array bounds and points, and x_0 x_1 x_2 on [0, 1]^3.
   auto box = [](const Real& c, const Real& a, const Real& b) { return (exp(c * b) - exp(c * a)) / c; };
   const std::array<Real, 3> a{0, -1, 3}, b{1, 2, 5};
   auto g = [](const std::array<Real, 3>& x) { return exp(x[0] / 2 - x[1] + x[2] / 4); };
   const Real G = sparse_grid(g, rule, a, b, 100000, tolerance, &error, &L1);
   check(box(Real(0.5), 0, 1) * box(Real(-1), -1, 2) * box(Real(0.25), 3, 5), G);
   // The result, error and L1 of the grid on [-1, 1]^3 for g(x(t)), times the volume / 8.
   const Real Gt = sparse_grid([&](const std::array<Real, 3>& t) {
      return g(std::array<Real, 3>{Real(0.5) + t[0] / 2, Real(0.5) + Real(1.5) * t[1], 4 + t[2]});
   }, rule, 3, 100000, tolerance, &error2, &L2);
   CHECK_ULP_CLOSE(Real(0.75) * Gt, G, 8);
   CHECK_ULP_CLOSE(Real(0.75) * error2, error, 8);
   CHECK_ULP_CLOSE(Real(0.75) * L2, L1, 8);
   check(Real(1) / 8, sparse_grid([](const std::vector<Real>& x) { return x[0] * x[1] * x[2]; }, rule, std::vector<Real>(3, 0),
                                  std::vector<Real>(3, 1), 1000, tolerance, &error, &L1));
   // A complex f and a reversed coordinate: the integrals of e^{ix} over [0, pi] and over [0, -pi/2] are 2i and -1 + i.
   const Complex S = sparse_grid([](const std::vector<Real>& x) { return Complex(cos(x[0] + x[1]), sin(x[0] + x[1])); }, rule,
                                 std::vector<Real>{0, 0}, std::vector<Real>{pi, -pi / 2}, 10000, tolerance, &error, &L1);
   CHECK_LE(Real(abs(S - Complex(0, 2) * Complex(-1, 1))), (std::max)(error, 64 * eps * L1));
   // A box of zero volume, even where the Jacobi exponents are negative.
   CHECK_EQUAL(sparse_grid(g, leja_jacobi<Real>(Real(-0.25), Real(-0.25)), a, std::array<Real, 3>{1, -1, 5}, 1000, tolerance, &error, &L1), Real(0));
   CHECK_EQUAL(error, Real(0));
   CHECK_EQUAL(L1, Real(0));
   const std::vector<Real> Zg = sparse_grid([](const std::vector<Real>& x) { return std::vector<Real>{x[0], x[1]}; }, rule,
                                            std::vector<Real>{0, 1}, std::vector<Real>{1, 1}, 1000, tolerance, &error, &L1);
   CHECK_EQUAL(Zg.size(), std::size_t(2));
   CHECK_EQUAL(Zg[0], Real(0));
   CHECK_EQUAL(Zg[1], Real(0));
   CHECK_EQUAL(error, Real(0));
   CHECK_EQUAL(L1, Real(0));
   CHECK_THROW(sparse_grid(g, rule, a, std::array<Real, 3>{1, 2, inf}, 1000), std::domain_error);
   CHECK_THROW(sparse_grid(g, rule, std::vector<Real>(3, 0), std::vector<Real>(2, 1), 1000), std::domain_error);
   CHECK_THROW(sparse_grid(g, rule, std::vector<Real>(2, 0), std::vector<Real>(2, 1), 1000), std::domain_error);
   CHECK_THROW(sparse_grid(g, rule, std::vector<Real>(2, 0), std::vector<Real>(2, 0), 1000), std::domain_error);
}

// Threads that share a rule, before it has any points, get what a rule of their own gives them.
void test_threads()
{
   const leja_hermite<double> shared;
   std::vector<double> Q(8);
   std::vector<std::thread> threads;
   for (unsigned t = 0; t < Q.size(); ++t)
      threads.emplace_back([&, t] { Q[t] = shared.integrate([t](double x) { return std::cos((t + 1) * x / 2); }, 1e-12); });
   for (std::thread& thread : threads)
      thread.join();
   // The compiler may contract the sums into fused multiply-adds differently where they are inlined.
   const leja_hermite<double> alone;
   for (unsigned t = 0; t < Q.size(); ++t)
      CHECK_ULP_CLOSE(Q[t], alone.integrate([t](double x) { return std::cos((t + 1) * x / 2); }, 1e-12), 2);
   const unsigned n = static_cast<unsigned>(shared.abscissa().size());
   for (unsigned i = 0; i < n; ++i)
      CHECK_EQUAL(shared.weights(n)[i], alone.weights(n)[i]);
}

void test_errors()
{
   CHECK_THROW(leja_laguerre<double>(-1.0), std::domain_error);
   CHECK_THROW(leja_jacobi<double>(0.5, -1.5), std::domain_error);
   CHECK_THROW(leja_jacobi<double>(0.0, std::numeric_limits<double>::quiet_NaN()), std::domain_error);
   const double inf = std::numeric_limits<double>::infinity();
   CHECK_THROW(leja_laguerre<double>(+inf), std::domain_error);
   CHECK_THROW(leja_jacobi<double>(inf, 0.0), std::domain_error);
   CHECK_THROW(leja_jacobi<double>(0.0, inf), std::domain_error);
   leja_legendre<double> rule;
   CHECK_THROW(rule.weights(0), std::domain_error);
   CHECK_THROW(rule.weights(rule.max_points() + 1), std::domain_error);
   // A NaN at the fifth point, 0.66, reaches the error estimate.
   double error;
   rule.integrate([](double x) { return x > 0 && x < 1 ? std::numeric_limits<double>::quiet_NaN() : std::exp(x); }, 1e-10, &error);
   CHECK_NAN(error);
   auto one = [](const std::vector<double>&) { return 1.0; };
   CHECK_THROW(boost::math::quadrature::sparse_grid(one, rule, 0, 10), std::domain_error);
   CHECK_THROW(boost::math::quadrature::sparse_grid(one, rule, 2, 0), std::domain_error);
   CHECK_THROW(boost::math::quadrature::sparse_grid([](const std::array<double, 3>&) { return 1.0; }, rule, 2, 10), std::domain_error);
}

// With domain errors ignored, an invalid n or exponent gives no weights and a NaN integral.
void test_ignore_error()
{
   using ignore = boost::math::policies::policy<boost::math::policies::domain_error<boost::math::policies::ignore_error>>;
   auto x = [](double t) { return t; };
   double error, L1;
   const leja_legendre<double, ignore> rule;
   CHECK_EQUAL(rule.weights(0).size(), std::size_t(0));
   CHECK_EQUAL(rule.weights(rule.max_points() + 1).size(), std::size_t(0));
   const leja_laguerre<double, ignore> laguerre(-1.0);
   const leja_jacobi<double, ignore> jacobi(0.5, -1.5);
   CHECK_EQUAL(laguerre.weights(1).size(), std::size_t(0));
   CHECK_EQUAL((leja_jacobi<double, ignore>(std::numeric_limits<double>::infinity(), 0.0).weights(1).size()), std::size_t(0));
   CHECK_NAN(laguerre.integrate(x, 1e-10, nullptr, &L1));
   CHECK_NAN(L1);
   CHECK_NAN(jacobi.integrate(x, 1e-10, &error, &L1));
   CHECK_NAN(error);
   CHECK_NAN(L1);
   CHECK_EQUAL(jacobi.abscissa().size(), std::size_t(0));
}

int main()
{
   test_exactness<float>();
   test_exactness<double>();
   test_exactness<long double>();
   test_exactness<boost::multiprecision::cpp_bin_float_50>(20);

   test_nested<float>();
   test_nested<double>();

   test_same_points<float, leja_legendre>(128);
   test_same_points<float, leja_hermite>(100);
   test_same_points<long double, leja_legendre>(128);
   test_same_points<long double, leja_hermite>(100);
   test_same_points<boost::multiprecision::cpp_bin_float_50, leja_legendre>(40);
   test_same_points<boost::multiprecision::cpp_bin_float_50, leja_hermite>(40);
   test_same_points<float, leja_laguerre>(100);
   test_same_points<long double, leja_laguerre>(100);
   test_same_points<boost::multiprecision::cpp_bin_float_50, leja_laguerre>(40);
   test_same_points<float, leja_laguerre>(100, 2.5);
   test_same_points<long double, leja_jacobi>(100, -0.5, 2.0);
   test_same_points<boost::multiprecision::cpp_bin_float_50, leja_jacobi>(40, 0.5, -0.75);

   test_smooth<double>(30);
   // long double is quad precision on some platforms, aarch64 Linux among them, and needs more points.
   const bool quad = std::numeric_limits<long double>::digits > 64;
   test_smooth<long double>(quad ? 40 : 30);

   test_adaptive<double>();
   test_adaptive<long double>();

   test_complex<float>(20, 30);
   test_complex<double>(20, 30);
   test_complex<long double>(quad ? 30 : 20, quad ? 50 : 30);
   test_complex<boost::multiprecision::cpp_bin_float_50, boost::multiprecision::cpp_complex_50>(36, 64);

   test_vector<float>();
   test_vector<double>();
   test_vector<boost::multiprecision::cpp_bin_float_50, boost::multiprecision::cpp_complex_50>();

   test_sparse_grid<float>(1e-4f);
   test_sparse_grid<double>(1e-10);
   test_sparse_grid<boost::multiprecision::cpp_bin_float_50, boost::multiprecision::cpp_complex_50>(1e-12);
   test_vanishing();

   test_bounds<float>(1e-4f);
   test_bounds<double>(1e-10);
#ifdef BOOST_HAS_FLOAT128
   test_bounds<boost::multiprecision::float128, boost::multiprecision::complex128>(1e-20);
#endif
   test_bounds<boost::multiprecision::cpp_bin_float_50, boost::multiprecision::cpp_complex_50>(1e-20);

   test_threads();
   test_errors();
   test_ignore_error();
   return boost::math::test::report_errors();
}
