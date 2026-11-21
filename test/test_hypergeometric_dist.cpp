// Copyright John Maddock 2008
// Copyright Paul A. Bristow 2021
// Copyright Gautam Sewani 2008

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY throw_on_error
#include <boost/math/concepts/real_concept.hpp> // for real_concept
#include <boost/math/distributions/hypergeometric.hpp>

#include "math_unit_test.hpp"

#include <iostream>
#include <cstdlib>
   using std::cout;
   using std::endl;
   using std::setprecision;

#include <array>
#include <boost/math/special_functions/relative_difference.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/stats.hpp>
#include <iomanip>
#include <stdexcept>
#include "table_type.hpp"
#include <cstdint>

#ifdef _MSC_VER
#pragma warning (disable:4127 4512)
#elif __GNUC__ >= 5
#  pragma GCC diagnostic ignored "-Woverflow"
#elif defined(__clang__)
#  pragma clang diagnostic ignored "-Wliteral-range"
#endif

#define CHECK_EX(a) \
   {\
      if(!CHECK_TRUE(a))\
      {\
         std::cerr << "Failure was with data ";\
         std::cerr << std::setprecision(35); \
         std::cerr << "x = " << x << ", r = " << r << ", n = " << n\
         << ", N = " << N << ", p = " << cp << ", q = " << ccp << std::endl;\
      }\
   }

struct error_bound { double peak, mean; };

// The largest errors expected of each type, in units of epsilon:
template <class Real>
error_bound expected_error(bool random_data)
{
   const bool is_concept = std::is_same<Real, boost::math::concepts::real_concept>::value;
   const bool is_double = std::is_same<Real, double>::value;
   bool is_largest = is_concept || std::is_same<Real, long double>::value;
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   using namespace boost::math::policies;
   if(!promote_double<>() || digits<double, policy<> >() == digits<long double, policy<> >())
   {
      is_largest = is_largest || is_double;
   }
   else if(is_double && (boost::math::tools::digits<long double>() < 100))
   {
      // Some split of errors from long double into double:
      return random_data ? error_bound{ 1500, 1500 } : error_bound{ 10, 10 };
   }
#else
   is_largest = is_largest || is_double;
#endif
   if(random_data && is_concept)
   {
      return { 250000000, 25000000 };
   }
   if(is_largest)
   {
      return random_data ? error_bound{ 10000000, 5000000 } : error_bound{ 50, 20 };
   }
   return { 1, 1 };
}

template <class T>
inline std::uint64_t make_unsigned(T x)
{
   return static_cast<std::uint64_t>(x);
}
template<>
inline std::uint64_t make_unsigned(boost::math::concepts::real_concept x)
{
   return static_cast<std::uint64_t>(x.value());
}

template <class T>
T pdf_tester(T r, T n, T N, T x)
{
   boost::math::hypergeometric_distribution<T> d(make_unsigned(r), make_unsigned(n), make_unsigned(N));
   return pdf(d, x);
}

template <class T>
T cdf_tester(T r, T n, T N, T x)
{
   boost::math::hypergeometric_distribution<T> d(make_unsigned(r), make_unsigned(n), make_unsigned(N));
   return cdf(d, x);
}

template <class T>
T ccdf_tester(T r, T n, T N, T x)
{
   boost::math::hypergeometric_distribution<T> d(make_unsigned(r), make_unsigned(n), make_unsigned(N));
   return cdf(complement(d, x));
}

// Checks f against column result_column of the table, whose first four columns are its arguments:
template <class Real, class T, class F>
void test_table(const T& data, F f, unsigned result_column, const char* type_name, const char* function_name, error_bound bound)
{
   boost::math::tools::stats<Real> errors;
   std::size_t worst = 0;
   for(std::size_t i = 0; i < data.size(); ++i)
   {
      const Real expected = data[i][result_column];
      Real computed = expected;
      try
      {
         computed = f(Real(data[i][0]), Real(data[i][1]), Real(data[i][2]), Real(data[i][3]));
      }
      catch(const std::underflow_error&)
      {
         computed = 0;
      }
      catch(const std::overflow_error&)
      {
         computed = boost::math::tools::max_value<Real>();
      }
      catch(const std::exception& e)
      {
         CHECK_TRUE(false);
         std::cerr << "  " << function_name << "<" << type_name << "> row " << i << ": " << e.what() << "\n";
      }
      const Real error = boost::math::relative_difference(computed, expected);
      if(!(error <= (errors.max)()))
      {
         worst = i;
      }
      // A non-finite result where a finite one is expected, or a gross error, fails by itself:
      if(!CHECK_TRUE(error <= 0.5f))
      {
         std::cerr << "  " << function_name << "<" << type_name << "> row " << i << ": found " << std::setprecision(35) << computed << ", expected " << expected << "\n";
      }
      errors.add(error);
   }
   const Real eps = boost::math::tools::epsilon<Real>();
   const Real peak = (errors.max)() / eps;
   const Real mean = errors.rms() / eps;
   std::cout << std::setprecision(4) << function_name << "<" << type_name << "> Max = " << peak << " RMS Mean=" << mean << ", worst at row " << worst << std::endl;
   CHECK_TRUE(peak <= bound.peak);
   CHECK_TRUE(mean <= bound.mean);
}

template <class Real, class T>
void do_test_hypergeometric(const T& data, const char* type_name, const char* test_name, error_bound bound)
{
   std::cout << "Testing " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   test_table<Real>(data, pdf_tester<Real>, 4, type_name, "hypergeometric PDF", bound);
   test_table<Real>(data, cdf_tester<Real>, 5, type_name, "hypergeometric CDF", bound);
   test_table<Real>(data, ccdf_tester<Real>, 6, type_name, "hypergeometric CDF complement", bound);
   std::cout << std::endl;
}

template <class Real, class T>
void do_test_hypergeometric_quantile(const T& data, const char* type_name, const char* test_name)
{
   typedef Real                   value_type;

   std::cout << "Checking quantiles with " << test_name << " with type " << type_name
      << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

   if(boost::math::tools::digits<value_type>() > 50)
   {
      for(unsigned i = 0; i < data.size(); ++i)
      {
         using namespace boost::math::policies;

         unsigned r = make_unsigned(data[i][0]);
         unsigned n = make_unsigned(data[i][1]);
         unsigned N = make_unsigned(data[i][2]);
         unsigned x = make_unsigned(data[i][3]);
         value_type cp = data[i][5];
         value_type ccp = data[i][6];

         boost::math::hypergeometric_distribution<value_type,
            policy<discrete_quantile<integer_round_up> > > du(r, n, N);

         if((cp < 0.9) && (cp > boost::math::tools::min_value<value_type>()))
         {
            CHECK_EX(quantile(du, cp) >= x);
         }
         if((ccp < 0.9) && (ccp > boost::math::tools::min_value<value_type>()))
         {
            CHECK_EX(quantile(complement(du, ccp)) >= x);
         }
         boost::math::hypergeometric_distribution<value_type,
            policy<discrete_quantile<integer_round_down> > > dl(r, n, N);
         if((cp < 0.9) && (cp > boost::math::tools::min_value<value_type>()))
         {
            CHECK_EX(quantile(dl, cp) <= x);
         }
         if((ccp < 0.9) && (ccp > boost::math::tools::min_value<value_type>()))
         {
            CHECK_EX(quantile(complement(dl, ccp)) <= x);
         }
         boost::math::hypergeometric_distribution<value_type,
            policy<discrete_quantile<integer_round_nearest> > > dn(r, n, N);

         if((cp < 0.9) && (cp > boost::math::tools::min_value<value_type>()))
         {
            CHECK_EX(quantile(dn, cp) == x);
         }
         if((ccp < 0.9) && (ccp > boost::math::tools::min_value<value_type>()))
         {
            CHECK_EX(quantile(complement(dn, ccp)) == x);
         }
         boost::math::hypergeometric_distribution<value_type,
            policy<discrete_quantile<integer_round_outwards> > > dou(r, n, N);

         if((cp < 0.9) && (cp > boost::math::tools::min_value<value_type>()))
         {
            if(cp < 0.5)
            {
               CHECK_EX(quantile(dou, cp) <= x);
            }
            else
            {
               CHECK_EX(quantile(dou, cp) >= x);
            }
         }
         if((ccp < 0.9) && (ccp > boost::math::tools::min_value<value_type>()))
         {
            if(ccp < 0.5)
            {
               CHECK_EX(quantile(complement(dou, ccp)) >= x);
            }
            else
            {
               CHECK_EX(quantile(complement(dou, ccp)) <= x);
            }
         }
         boost::math::hypergeometric_distribution<value_type,
            policy<discrete_quantile<integer_round_inwards> > > di(r, n, N);

         if((cp < 0.9) && (cp > boost::math::tools::min_value<value_type>()))
         {
            if(cp < 0.5)
            {
               CHECK_EX(quantile(di, cp) >= x);
            }
            else
            {
               CHECK_EX(quantile(di, cp) <= x);
            }
         }
         if((ccp < 0.9) && (ccp > boost::math::tools::min_value<value_type>()))
         {
            if(ccp < 0.5)
            {
               CHECK_EX(quantile(complement(di, ccp)) <= x);
            }
            else
            {
               CHECK_EX(quantile(complement(di, ccp)) >= x);
            }
         }
      }
   }
}


template <class RealType>
void test_spot(unsigned x, unsigned n, unsigned r, unsigned N,
               RealType p, RealType cp, RealType ccp, RealType tol)
{
   boost::math::hypergeometric_distribution<RealType> d(r, n, N);

   std::pair<unsigned, unsigned> extent = range(d);
   // CDF's:
   CHECK_RELATIVE_ERROR(pdf(d, x), p, tol);
   CHECK_RELATIVE_ERROR(cdf(d, x), cp, tol);
   CHECK_RELATIVE_ERROR(cdf(complement(d, x)), ccp, tol);
   // Again with real-value arguments:
   CHECK_RELATIVE_ERROR(pdf(d, static_cast<RealType>(x)), p, tol);
   CHECK_RELATIVE_ERROR(cdf(d, static_cast<RealType>(x)), cp, tol);
   CHECK_RELATIVE_ERROR(cdf(complement(d, static_cast<RealType>(x))), ccp, tol);
   //
   // Quantiles, don't bother checking these for type float
   // as there's not enough precision in the p and q values
   // to get back to where we started:
   //
   if(boost::math::tools::digits<RealType>() > 50)
   {
      using namespace boost::math::policies;

      boost::math::hypergeometric_distribution<RealType,
         policy<discrete_quantile<integer_round_up> > > du(r, n, N);

      CHECK_EX(quantile(du, cp) >= x);
      CHECK_EX(quantile(complement(du, ccp)) >= x);

      boost::math::hypergeometric_distribution<RealType,
         policy<discrete_quantile<integer_round_down> > > dl(r, n, N);
      CHECK_EX(quantile(dl, cp) <= x);
      CHECK_EX(quantile(complement(dl, ccp)) <= x);

      boost::math::hypergeometric_distribution<RealType,
         policy<discrete_quantile<integer_round_nearest> > > dn(r, n, N);

      CHECK_EX(quantile(dn, cp) == x);
      CHECK_EX(quantile(complement(dn, ccp)) == x);
   }
   //
   // Error checking of out of bounds arguments:
   //
   CHECK_THROW(pdf(d, extent.second + 1), std::domain_error);
   CHECK_THROW(cdf(d, extent.second + 1), std::domain_error);
   CHECK_THROW(cdf(complement(d, extent.second + 1)), std::domain_error);
   if(extent.first > 0)
   {
      CHECK_THROW(pdf(d, extent.first - 1), std::domain_error);
      CHECK_THROW(cdf(d, extent.first - 1), std::domain_error);
      CHECK_THROW(cdf(complement(d, extent.first - 1)), std::domain_error);
   }
   CHECK_THROW(quantile(d, 1.1f), std::domain_error);
   CHECK_THROW(quantile(complement(d, 1.1f)), std::domain_error);
   CHECK_THROW(quantile(d, -0.001f), std::domain_error);
   CHECK_THROW(quantile(complement(d, -0.001f)), std::domain_error);
   //
   // Checking of extreme values:
   //
   CHECK_EQUAL(quantile(d, 0), RealType(extent.first));
   CHECK_EQUAL(quantile(d, 1), RealType(extent.second));
   CHECK_EQUAL(quantile(complement(d, 0)), RealType(extent.second));
   CHECK_EQUAL(quantile(complement(d, 1)), RealType(extent.first));
   CHECK_EQUAL(cdf(d, extent.second), RealType(1));
   CHECK_EQUAL(cdf(complement(d, extent.second)), RealType(0));
}

template <class RealType>
void test_spots(RealType /*T*/, const char* type_name)
{
   // Basic sanity checks.
   // 50 eps as a percentage, up to a maximum of double precision
   // Test data taken from Mathematica 6
// The tables' SC_ uses BOOST_JOIN, from Boost.Config, which a standalone build doesn't have:
#ifndef BOOST_JOIN
#define BOOST_JOIN(X, Y) BOOST_MATH_JOIN(X, Y)
#endif
#define T RealType
#include "hypergeometric_test_data.ipp"
   do_test_hypergeometric<T>(hypergeometric_test_data, type_name, "Mathematica data", expected_error<T>(false));

#include "hypergeometric_dist_data2.ipp"
   if(std::is_floating_point<RealType>::value)
   {
      //
      // Don't test this for real_concept: it's too slow!!!
      //
      do_test_hypergeometric<T>(hypergeometric_dist_data2, type_name, "Random large data", expected_error<T>(true));
   }

   do_test_hypergeometric_quantile<T>(hypergeometric_test_data, type_name, "Mathematica data");
   if(std::is_floating_point<RealType>::value)
   {
      //
      // Don't test this for real_concept: it's too slow!!!
      //
      do_test_hypergeometric_quantile<T>(hypergeometric_dist_data2, type_name, "Random large data");
   }

   RealType tolerance = (std::max)(
      static_cast<RealType>(2e-16L),    // limit of test data
      boost::math::tools::epsilon<RealType>());
   cout<<"Absolute tolerance:"<<tolerance<<endl;

   tolerance *= 50;
   cout << "Tolerance for type " << typeid(RealType).name()  << " is " << tolerance << endl;

   //
   // These sanity check values were calculated using the online
   // calculator at http://stattrek.com/Tables/Hypergeometric.aspx
   // It's assumed that the test values are accurate to no more than
   // double precision.
   //
   test_spot(20, 200, 50, 500, static_cast<T>(0.120748236361163), static_cast<T>(0.563532430195156), static_cast<T>(1 - 0.563532430195156), tolerance);
   test_spot(53, 452, 64, 500, static_cast<T>(0.0184749573044286), static_cast<T>(0.0299118078796907), static_cast<T>(1 - 0.0299118078796907), tolerance);
   test_spot(32, 1287, 128, 5000, static_cast<T>(0.0807012167418264), static_cast<T>(0.469768774237742), static_cast<T>(1 - 0.469768774237742), tolerance);
   test_spot(1, 13, 4, 26, static_cast<T>(0.248695652173913), static_cast<T>(0.296521739130435), static_cast<T>(1 - 0.296521739130435), tolerance);
   test_spot(2, 13, 4, 26, static_cast<T>(0.40695652173913), static_cast<T>(0.703478260869565), static_cast<T>(1 - 0.703478260869565), tolerance);
   test_spot(3, 13, 4, 26, static_cast<T>(0.248695652173913), static_cast<T>(0.952173913043478), static_cast<T>(1 - 0.952173913043478), tolerance);
   test_spot(40, 70, 89, 170, static_cast<T>(0.0721901023798991), static_cast<T>(0.885447799131944), static_cast<T>(1 - 0.885447799131944), tolerance);

   boost::math::hypergeometric_distribution<RealType> d(50, 200, 500);
   CHECK_EQUAL(range(d).first, UINT64_C(0));
   CHECK_EQUAL(range(d).second, UINT64_C(50));
   CHECK_RELATIVE_ERROR(mean(d), static_cast<RealType>(20), tolerance);
   CHECK_RELATIVE_ERROR(mode(d), static_cast<RealType>(20), tolerance);
   // Test values and code revised to correct kurtosis using Mathematica algorithm and test values.
   // Also checked using boost::multiprecision::cpp_bin_float_quad at 128-bit precision.
   // And compared to an implementation of the Wikipedia equation at https://en.wikipedia.org/wiki/Hypergeometric_distribution.
   // https://github.com/boostorg/math/issues/639  hypergeometric kurtosis differs from Wolfram and scipy and Wikipedia.
   // N[variance[hypergeometricdistribution(200,50,500)], 55]
   CHECK_RELATIVE_ERROR(variance(d), static_cast<RealType>(10.82164328657314629258517034068136272545090180360721443L), tolerance);
   // N[skewness[hypergeometricdistribution(200,50,500)], 55]
   CHECK_RELATIVE_ERROR(skewness(d), static_cast<RealType>(0.04883307102295208473290291018936563084537649911337978351L), tolerance);
   // https://www.wolframalpha.com/input/?i=N%5Bkurtosis%5Bhypergeometricdistribution%28200%2C50%2C500%29%5D%2C+55%5D+
   // N[kurtosis[hypergeometricdistribution(200,50,500)], 55] = 2.969174035736058474901169623721804275002985337280263464
   CHECK_RELATIVE_ERROR(kurtosis(d), static_cast<RealType>(2.969174035736058474901169623721804275002985337280263464L), tolerance);
   CHECK_RELATIVE_ERROR(kurtosis(d), kurtosis_excess(d) + 3, tolerance);
   CHECK_EQUAL(quantile(d, 0.5f), median(d));

   CHECK_THROW(d = boost::math::hypergeometric_distribution<RealType>(501, 40, 500), std::domain_error);
   CHECK_THROW(d = boost::math::hypergeometric_distribution<RealType>(40, 501, 500), std::domain_error);

   //
   // Integer x beyond the factorial table, where Loader's saddle-point method is used, with prime
   // factorisation in the tails up to N = 104729. Values from mpmath at 60 digits.
   //
   struct large_case { unsigned x, r, n, N; long double pdf; RealType tol; };
   const large_case large_cases[] = {
      { 600, 1500, 2000, 5000, 0.025126012323870050456659927618198748L, tolerance },
      { 700, 1500, 2000, 5000, 0.0000000000708147414413862429733772355814014052L, tolerance },
      { 24000, 60000, 80000, 200000, 0.0039735405018628404013758812090728637L, tolerance },
      { 24123, 60000, 80000, 200000, 0.00187568978055549884111723163694289891L, tolerance },
      { 23500, 60000, 80000, 200000, 0.000000016112756785169958209568753951147937L, 4 * tolerance },
      { 1200000, 3000000, 4000000, 10000000, 0.000561946220334811372962916024931211627L, tolerance },
      { 1203456, 3000000, 4000000, 10000000, 0.00000000402250531175235295684236178111301572L, 4 * tolerance },
   };
   for (auto const& c : large_cases)
   {
      boost::math::hypergeometric_distribution<RealType> dist(c.r, c.n, c.N);
      CHECK_RELATIVE_ERROR(pdf(dist, c.x), static_cast<RealType>(c.pdf), c.tol);
   }
   // At an end of the support the result must not underflow before it is complete (true value 1.466e-322).
   if (std::numeric_limits<RealType>::min_exponent10 < -320)
   {
      boost::math::hypergeometric_distribution<RealType> dist(100000, 1065, 200000);
      CHECK_TRUE(pdf(dist, 0u) > 0);
   }

   //
   // Non-integer x: the analytic continuation C(r, x) C(N - r, n - x) / C(N, n).
   // Values from mpmath at 50 digits; the tolerance scales with the condition number
   // |x d/dx ln pdf|, which is how much a rounding of x itself is amplified.
   //
   struct noninteger_case { const char* x; unsigned r, n, N; long double pdf, cond; };
   static const noninteger_case noninteger_cases[] = {
      { "2.5", 7, 9, 20, 0.306765408837835280425772348087512074L, 1.27L },
      { "0.25", 7, 9, 20, 0.00923344250444341451195317576706156765L, 0.72L },
      { "6.75", 7, 9, 20, 0.00121681025252355074963724021736470271L, 24.4L },
      { "60.5", 200, 150, 500, 0.0788614625136665537116488597232431263L, 1.29L },
      { "45.125", 200, 150, 500, 0.000938770338227857359230205729668782774L, 27.4L },
      { "12000.5", 30000, 40000, 100000, 0.00561924531383347907377494860263026125L, 1.29L },
   };
   for (auto const& c : noninteger_cases)
   {
      boost::math::hypergeometric_distribution<RealType> dist(c.r, c.n, c.N);
      const RealType x = static_cast<RealType>(std::strtold(c.x, nullptr));
      CHECK_RELATIVE_ERROR(pdf(dist, x), static_cast<RealType>(c.pdf), tolerance * static_cast<RealType>(1 + c.cond));
   }
   // Integer-valued x of any type takes the exact integer path.
   boost::math::hypergeometric_distribution<RealType> d2(200, 150, 500);
   for (unsigned k = 0; k <= 150; k += 10)
   {
      CHECK_EQUAL(pdf(d2, static_cast<RealType>(k)), pdf(d2, static_cast<std::uint64_t>(k)));
   }
   // The continuation is log-concave in x, so between integers away from the mode it is monotone.
   CHECK_TRUE(pdf(d2, static_cast<RealType>(40.5)) > pdf(d2, static_cast<RealType>(40)));
   CHECK_TRUE(pdf(d2, static_cast<RealType>(40.5)) < pdf(d2, static_cast<RealType>(41)));
   CHECK_TRUE(pdf(d2, static_cast<RealType>(80.5)) < pdf(d2, static_cast<RealType>(80)));
   CHECK_TRUE(pdf(d2, static_cast<RealType>(80.5)) > pdf(d2, static_cast<RealType>(81)));
   // Counts too large for the type to hold exactly (here N > 2^24 for an unpromoted float): their
   // differences must be formed exactly, or at the bottom of the support (7) a gamma argument hits a pole.
   {
      using no_promote = boost::math::policies::policy<boost::math::policies::promote_float<false>, boost::math::policies::promote_double<false>>;
      boost::math::hypergeometric_distribution<RealType, no_promote> big(30000001 - 3, 10, 30000001);
      const RealType p = pdf(big, static_cast<RealType>(7.5));
      CHECK_TRUE((boost::math::isfinite)(p) && (p >= 0));
   }
   // Outside the support is still an error.
   CHECK_THROW(pdf(d2, static_cast<RealType>(150.5)), std::domain_error);
   CHECK_THROW(pdf(d2, static_cast<RealType>(-0.5)), std::domain_error);
}


template <class RealType>
void test_entropy(RealType)
{
   // Reference values from mpmath at 60 digits: closed forms, or -sum p log(p) for the discrete
   // distributions, which have none. Parameters are exact in every type.
   RealType tol = 16 * boost::math::tools::epsilon<RealType>();
   auto check = [&](std::uint64_t a, std::uint64_t b, std::uint64_t c, RealType h)
   {
      CHECK_RELATIVE_ERROR(boost::math::entropy(boost::math::hypergeometric_distribution<RealType>(a, b, c)), h, tol);
   };
   check(30, 40, 100, static_cast<RealType>(2.23251995010132138603305010352835455335946637L));
   check(5, 10, 20, static_cast<RealType>(1.41017899258903122865281944730670874143759777L));
   check(1, 1, 2, static_cast<RealType>(0.693147180559945309417232121458176568075500134L));
   check(10, 1000, 1000000, static_cast<RealType>(0.0560825642231663484256230658066717444208256013L));
   check(300000, 400000, 1000000, static_cast<RealType>(6.83281224965225896648567487522138918484639344L));
   check(250000000, 333333333, 1000000000, static_cast<RealType>(10.1915445366983499530007834277413067488348034L));
   CHECK_EQUAL(boost::math::entropy(boost::math::hypergeometric_distribution<RealType>(5, 5, 5)), RealType(0));
}

int main()
{
   // Basic sanity-check spot values.
   // (Parameter value, arbitrarily zero, only communicates the floating point type).
   test_spots(0.0F, "float"); // Test float. OK at decdigits = 0 tolerance = 0.0001 %
   test_entropy(0.0F);
   test_spots(0.0, "double"); // Test double. OK at decdigits 7, tolerance = 1e07 %
   test_entropy(0.0);
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   test_spots(0.0L, "long double"); // Test long double.
   test_entropy(0.0L);
#ifndef BOOST_MATH_NO_REAL_CONCEPT_TESTS
   test_spots(boost::math::concepts::real_concept(0), "real_concept"); // Test real_concept.
   test_entropy(boost::math::concepts::real_concept(0));
#endif
#else
   std::cout << "<note>The long double tests have been disabled on this platform "
      "either because the long double overloads of the usual math functions are "
      "not available at all, or because they are too inaccurate for these tests "
      "to pass.</note>" << std::endl;
#endif

   return boost::math::test::report_errors();
}

