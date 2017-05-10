// Copyright Thomas Luu, 2014
// Copyright Paul A. Bristow 2016, 2017

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or
//  copy at http ://www.boost.org/LICENSE_1_0.txt).

// Computation of the Lambert W function using the algorithm from
// Thomas Luu, UCL Department of Mathematics, PhD Thesis, (2016)
// Fast and Accurate Parallel Computation of Quantile Functions for Random Number Generation,
// page 96 - 98, text and Routine 11.

// This implementation of the algorithm computes W(x) when x is any C++ real type,
// including Boost.Multiprecision types, and also the prototype Boost.Fixed_point types.

// This version only handles real values of W(x) when x > -1/e (x > -0.367879...),
// and just returns NaN for x < -1/e.

// Initial guesses based on :
// D.A.Barry, J., Y.Parlange, L.Li, H.Prommer, C.J.Cunningham, and F.Stagnitti.
// Analytical approximations for real values of the Lambert W function.
// Mathematics and Computers in Simulation, 53(1) : 95 - 103, 2000.

// D.A.Barry, J., Y.Parlange, L.Li, H.Prommer, C.J.Cunningham, and F.Stagnitti.
// Erratum to analytical approximations for real values of the Lambert W function.
// Mathematics and computers in simulation, 59(6) : 543 - 543, 2002.

// See also https://svn.boost.org/trac/boost/ticket/11027
// and discussions at http://lists.boost.org/Archives/boost/2016/09/230832.php
// This code gives identical values as https://github.com/thomasluu/plog
// https://github.com/thomasluu/plog/blob/master/plog.cu
// for type double.


#ifndef BOOST_MATH_SF_LAMBERT_W_HPP
#define BOOST_MATH_SF_LAMBERT_W_HPP

//#define BOOST_MATH_INSTRUMENT_LAMBERT_W  // Only #define for Lambert W diagnostic output.

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/log1p.hpp> // for log (1 + x)
#include <boost/math/constants/constants.hpp> // For exp_minus_one == 3.67879441171442321595523770161460867e-01.
#include <boost/math/special_functions/next.hpp> // for float_distance
#include <boost/math/policies/policy.hpp>

#include <boost/fixed_point/fixed_point.hpp> // fixed_point, if required.
#if defined (BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable: 4127) // Conditional expression is constant.
#endif


int total_iterations = 0;
int total_calls = 0;
int max_iterations = 0;

double worst_diff = 0;
// boost::math::tools::max_value<double>();
double total_diffs = 0.;

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
#  include <iostream>  // Only needed for testing.
#endif
#include <limits> // numeric_limits infinity & NaN, for diagnostic max_digits10 and digits10 ...
#include <cmath> // exp function.
#include <type_traits> // is_integral trait.

namespace boost {
  namespace fixed_point {
    //! \brief Specialization for log1p for fixed_point types (in its own fixed_point namespace).

    /*!
    \details This simple (and fast) approximation gives a result
    within a few epsilon/ULP of the nearest representable value.
    This makes it suitable for use as a first approximation by the Lambert W function below.
    It only uses (and thus relies on for accuracy of) the standard log function.
    Refinement of the Lambert W estimate using Halley's method
    seems to get quickly to high accuracy in a very few iterations,
    so that any small inaccuracy in the 1st approximation is fairly unimportant.
    It also avoids any Taylor series refinements that involve high powers of x
    that could easily go outside the often limited range of fixed-point types.
    */

    template<const int IntegralRange, const int FractionalResolution, typename RoundMode, typename OverflowMode>
    negatable<IntegralRange, FractionalResolution, RoundMode, OverflowMode>
      log1p(negatable<IntegralRange, FractionalResolution, RoundMode, OverflowMode> x)
    {
      // https://github.com/f32c/arduino/blob/master/hardware/fpga/f32c/system/src/math/log1p.c
      // HP - 15C Advanced Functions Handbook, p.193.

      typedef negatable<IntegralRange, FractionalResolution, RoundMode, OverflowMode> local_negatable_type;

      local_negatable_type unity(1);
      // A fixed_point type might not include unity: hope something sensible happens.
      local_negatable_type u = unity + x;
      if (u == unity)
      { //
        return x;
      }
      else
      {
        local_negatable_type result = log(u) * (x / (u - unity));
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
        std::cout << "log1p fixed-point approx " << result << std::endl;
#endif
        return result;
      } // if else
    } // log1p_fp

  } //namespace fixed_point
} // namespace boost

namespace local
{
  // Special versions of log1p that may be useful.

  // log1p_dodgy version seems to work OK, but...
  /*  John Maddock cautions that
  "The formula relies on:

  1) That the operations are not re-ordered.

  2) That each operation is correctly rounded.
  Is true only for floating-point arithmetic unless you compile with
  special flags to force full-IEEE compatible mode,
  but that slows the whole program down...

  3) That compiler optimisations don't perform symbolic math on the expressions and simplify them."
  */

  template<typename T>
  T log1p_dodgy(T x)
  {
    // log(1+x) == (log(1+x) * x) / ((1-x) - 1)
    T result = -(log(1 + x) * x) / ((1 - x) - 1); // From note in Boost.Math log1p.

    // \boost\libs\math\include\boost\math\special_functions\log1p.hpp
    // Algorithm log1p is part of C99, but is not yet provided by all compilers.
    //
    // The Boost.Math version uses a Taylor series expansion for 0.5 > x > epsilon, which may
    // require up to std::numeric_limits<T>::digits+1 terms to be calculated.
    // It would be much more efficient to use the equivalence:
    //   log(1+x) == (log(1+x) * x) / ((1-x) - 1)
    // Unfortunately many optimizing compilers make such a mess of this, that
    // it performs no better than log(1+x): which is to say not very well at all.

    return result;
  } // template<typename T> T log1p_dodgy(T x)

  template<typename T>
  inline
    T log1p(T x)
  {
    //! \brief log1p function.
    //! This is probably faster than using boost::math::log1p (if a little less accurate)
    // but should be good enough for computing Lambert_w initial estimate.

    //! \details Formula for function log1p from a Note in
    // https://github.com/f32c/arduino/blob/master/hardware/fpga/f32c/system/src/math/log1p.c
    // HP - 15C Advanced Functions Handbook, p.193.
    // Assuming log() returns an accurate answer,
    // the following algorithm can be used to compute log1p(x) to within a few ULP :
    // u = 1 + x;
    //  if (u == 1) return x
    // else return log(u)*(x/(u-1.));
    // This implementation is template of type T

    // log(u) * (x / (u - unity));

    // http://thepeg.hepforge.org/svn/tags/RELEASE_1_0/ThePEG/Utilities/Math.icc
    //   double log1m(double x) { return log1p(-x);}

    // It is probably possible to use Newton's or Halley's method to refine this approximation,
    // but improvement is probably not useful for estimating the Lambert W value because it is
    // close enough for the Halley's method to converge just as fast as with a better initial estimate.

    T unity;
    unity = static_cast<T>(1);
    T u = unity + x;
    if (u == unity)
    {
      return x;
    }
    else
    {
      T result = log(u) * (x / (u - unity));
      //T dodgy = log1p_dodgy(x);
      //double dx = static_cast<double>(x);
      //double lp1x = log1p(dx);
      /*
      #ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
      std::cout << "true " << lp1x
      << ", log1p_fp " << log1p(x)
      //<< ", dodgy " << dodgy << ", HP approximation " << result
      //<< ", diff " << dodgy - result << "epsilons = " << (dodgy - result)/std::numeric_limits<T>::epsilon()
      << std::endl;
      // For example: x = 0.5, dodgy 0.788421631, HP 0.78843689, diff -1.52587891e-05
      #endif
      */
      return result;
    }
  }  // template<typename T> T log1p(T x)
} // namespace local

namespace boost
{
  namespace math
  {
    namespace detail
    {
      //! \brief Approximation of real-branched values of the Lambert W function.
      //! \details Based on Barry et al, and Thomas Luu, Thesis (2016).
      //! Notes refer to Luu equations, page 96-98 and C++ code from algorithm Routine 11.
      // Assumes parameter x has been checked before call.
      // This function might be useful if the user only requires an approximate result quickly,
      // as relative errors of 0.02% or less are claimed by Barry et al.

      template <typename RealType = float>
      inline BOOST_CONSTEXPR_OR_CONST
      RealType lambert_w_approx(RealType x)
      {
        BOOST_MATH_STD_USING  // for ADL of std functions.

        using boost::math::constants::e; // e(1) = 2.71828...
        using boost::math::constants::root_two; // 1.414 ...
        using boost::math::constants::one_div_root_two; // 0.707106

        //using local::log1p; // Simplified approximation for speed. (mean 0.222 us)
        using boost::math::log1p; // Other approximate implementations may also be used.

        RealType w0 = 0; // TODO avoid this assignment to avoid GCC unitialized variable 'w0' in constexpr function.
        // Upper branch W0 is divided in two regions, W0+ and W0-: -1 <= w0_minus <= 0 and 0 <= w0_plus.
        if (x > 0)
        { // W0+ branch, Luu Routine 11, line 8, and equation 6.44, from Barry et al.
          // (Factors 1.2 and 2.4 are 'exact' from integral fractions 6/5 and 12/5).
          // A maximum relative error of 0.196% is claimed by Barry for equation 15, page 7.
          // Luu claims a maximum relative error of 2.39% for this slightly simplified version.
          w0 = log(static_cast<RealType>(1.2L) * x / log(static_cast<RealType>(2.4L) * x / log1p(static_cast<RealType>(2.4L) * x)));
          // local::log1p does seem to refine OK with multiprecision.
        }
        else
        { // Lower branch W0-> -exp(-1) or > -0.367879, so for a Real result need x > -0.367879 >= 0.
          // Luu, equation 6.44, page 97 claims a maximum relative error of 0.013%.
          // from Barry et al (Erratum Mathematics and Computers in Simulation, 53 (2000) 95–103)
          // Section 3.2 Approximation for W0-, equation 9 (but corrected in the erratum).
          // Barry et al claim a maximum relative error of 0.013%.
          RealType sqrt_v = root_two<RealType>() * sqrt(static_cast<RealType>(1) + e<RealType>() * x); // nu = sqrt(2 + 2e * z) Line 10.
          RealType n2 = static_cast<RealType>(3) * root_two<RealType>() + static_cast<RealType>(6); // 3 * sqrt(2) + 6, Line 11.
                                                                                                    // == 10.242640687119285146
          RealType n2num = (static_cast<RealType>(2237) + (static_cast<RealType>(1457) * root_two<RealType>())) * e<RealType>()
            - (static_cast<RealType>(4108) * root_two<RealType>()) - static_cast<RealType>(5764); //  Numerator Line 11.
          RealType n2denom = (static_cast<RealType>(215) + (static_cast<RealType>(199) * root_two<RealType>())) * e<RealType>()
            - (430 * root_two<RealType>()) - static_cast<RealType>(796); // Denominator Line 11.
          RealType nn2 = n2num / n2denom; // Line 11 full fraction.
          nn2 *= sqrt_v; // Line 11 last part.
          n2 -= nn2; //  Line 11 complete, equation 6.44,
          RealType n1 = (static_cast<RealType>(1) - one_div_root_two<RealType>()) * (n2 + root_two<RealType>()); // Line 12.
          w0 = -1 + sqrt_v * (n2 + sqrt_v) / (n2 + sqrt_v + n1 * sqrt_v); // W0 1st approximation, Line 13, equation 6.40.
        }
        return w0;
      } // template <typename RealType> RealType lambert_w_approx(RealType x)

      //! \brief Lambert W approximation that can be used when x is large.
      //! This approximation from Darko Veberic, 2.1 Recursion, page 4, equation 14, (2009), https://arxiv.org/pdf/1003.1628.pdf.
      //! x > std::numeric_limits<float>::max(); FLOAT_MAX
      //! and so using the lambert_w_approx<float> version for speed is no longer possible.
      //! x must be < 1/<RealType>epsilon
      //! For multiprecision types, RealType epsilon is so small (< 1e38) that
      //! FLOAT_MAX < 1/epsilon float: 3.4e+38 or 0x1.fffffep+127
      //! This is never true for fundamental types, including binary128 quad precision,
      //! epsilon = 2^-112 ~=1.93E-34, 1/epsilon = 5.19e33
      //! which is < 3.4e+38 and so using float is OK.

      template <typename T>
      inline BOOST_CONSTEXPR_OR_CONST
      T lambert_w_ln_approx(T z)
      {
        T l = log(z);
        return l - log(l); // 1st two terms of continued logarithm.
      } // template <typename T> T lambert_W_ln_approx(T z)


      //! \brief  Lambert W function implementation.
      //! \details Based on Thomas Luu, Thesis (2016).
      //! Notes refer to equations, page 96-98 and C++ code from algorithm Routine 11.
      template <typename RealType = double, class Policy = policies::policy<> >
      RealType lambert_w_imp(RealType x, const Policy& /* pol */)
      {
        BOOST_MATH_STD_USING  // for ADL of std functions.

        //using boost::math::constants::root_two; // 1.414 ...
        //using boost::math::constants::e; // e(1) = 2.71828...
        //using boost::math::constants::one_div_root_two; // 0.707106
        using boost::math::constants::exp_minus_one; // 0.36787944
        using boost::math::tools::max_value;

        // Catch the very common mistake of providing an integer value as parameter x to lambert_w, for example, lambert_w(1).
        // Need to ensure it is a floating-point type (of the desired type, float 1.F, double 1., or long double 1.L),
        // or static_casted integer, for example:  static_cast<float>(1) or static_cast<cpp_dec_float_50>(1).
        // Want to allow fixed_point types too, so do not just test for floating-point.
        BOOST_STATIC_ASSERT_MSG(!std::is_integral<RealType>::value, "Must be floating-point type (not integer type), for example: W(1.), not W(1)!");
        total_calls++; // Temporary for testing.

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
        std::cout << std::showpoint << std::endl; // Show all trailing zeros.
        std::cout.precision(std::numeric_limits<RealType>::max_digits10); // Show all possibly significant digits (17 for double).
                                                                          //std::cout.precision(std::numeric_limits<RealType>::digits10); // Show all significant digits (15 for double).
#endif
        // Check on range of x.
        if (x == 0)
        { // Special case of zero (for speed and to avoid log(0)and /0 etc).
          return static_cast<RealType>(0);
        }
        if (!boost::math::isfinite(x))
        {  // Check x is finite.
          return policies::raise_domain_error<RealType>(
            "boost::math::lambert_w<%1%>",
            "The parameter %1% is not finite, so the return value is NaN.",
            x,
            Policy());
        } //  (!boost::math::isfinite(x))

        // Check x is not so large that it would cause overflow,
        // and would throw exception "Error in function boost::math::log1p<double>(double): numeric overflow"
        if (x > boost::math::tools::max_value<RealType>() / 4) // Close to max == 4.4942328371557893e+307
        {
          return policies::raise_overflow_error<RealType>(
            "boost::math::lambert_w<%1%>",
            "The parameter %1% is too large, so the return value is NaN.",
            x, Policy());
          // Exception Error in function boost::math::lambert_w<long double>: The parameter 4.6180304784183997e+307 is too large, so the return value is NaN.
        }

        // Special case of singularity at -exp(-1)) // -0.3678794411714423215955237701614608674458111310
        // std::cout << "-exp(-1) = " << -expminusone<RealType>() << std::endl;
        // https://www.wolframalpha.com/input/?i=-exp(-1)&wal=header  N[-Exp[-1], 143]
        // -0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527
        // Added to Boost math constants as exp_minus_one
        // Can't use if (x < -exp(-1)) because possible 1-bit inaccuracy of exp means is inconsistent.
        if (x == -exp_minus_one<RealType>())
        {  // At Singularity -0.36787944117144233 == -0.36787944117144233 returned -1.0000000000000000
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
          std::cout << "At Singularity " << x << " == " << -exp_minus_one<RealType>() << " returned " << static_cast<RealType>(-1.) << std::endl;
#endif
          return static_cast<RealType>(-1);
        }
        else if (x < -exp_minus_one<RealType>()) // -0.3678794411714423215955237701614608674458111310L
        { // If x < -exp(-1)) then lambert_W(x) would be complex (not handled with this implementation).
          return policies::raise_domain_error<RealType>(
            "boost::math::lambert_w<%1%>",
            "The parameter %1% would give a complex result, so the return value is NaN.",
            x,
            // "Function boost::math::lambert_w<double>: The parameter -0.36787945032119751 would give a complex result, so the return value is NaN."
            Policy());
        }

        // Get a suitable starting approximation of lambert_w.
        // For speed, use float because the error of the approximation is more than float precision,
        // and the iteration will probably get to the result just as quickly.
        // For multiprecision types, this can almost half the execution time.
        // Some very precise multiprecision types (using more than about 150 bits)
        // will not be representable using float types,
        // so a continued log approximation will be used instead.

        using boost::math::detail::lambert_w_approx;
        using boost::math::policies::digits2;
        using boost::math::policies::digits10;
        using boost::math::policies::digits_base10;
        using boost::math::policies::get_epsilon;

        RealType epsilon = get_epsilon<RealType, Policy>();
        // default epsilon double 2.2204460492503131e-16
#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
        int d10 = policies::digits_base10<RealType, Policy>(); // policy template parameter digits10 in lambert_w(x, policy<digits10<1> >()
        int d2 = policies::digits<RealType, Policy>();
        std::cout << "digits10 " << d10 <<", digits_2 " << d2 << ", epsilon " << epsilon << std::endl;
#endif  // BOOST_MATH_INSTRUMENT_LAMBERT_W

        //int d2 = policies::digits<RealType, Policy>(); //
        //std::cout << "digits2 " << d2 << std::endl; // default 53 , for policy policy<digits10<3> > 13
        //int d10 = policies::digits_base10<RealType, Policy>();
        //std::cout << "digits10 " << d10 << std::endl; // 15 = 53 * 0.301
        //int d = policies::digits < RealType, policies::policy<> >();
        //std::cout << "digits <RealType, policies::policy<>() " << d << std::endl; // default = 53,  for policy policy<digits10<3> > = 3

        int iterations_required = 0;
        if (policies::digits<RealType, Policy>() != policies::digits<RealType, policies::policy<> >())
        { // User template parameter digits10 has specified digits != default (a reduced precision),
          // so just guestimate how many iterations should give us this precision:
          int initial_precision_bits = 6; // Should be how many bits correct in the the initial approximation w0.
          // If relative precision is about 2%, then about 2 decimal places, so very about 6 bits.
          // So should not do any iterations if < policies::digits<RealType, Policy>()
          // Can't be more than the binary precision of the RealType, say 53 for double.

          while (initial_precision_bits < policies::digits<RealType, Policy>())
          {
            initial_precision_bits *= 3;  // Assume each iteration doubles or triples the precision?
            // With Halley's method, tripling precision with each iteration might be more realistic?
            ++iterations_required;
          } // while
        //  std::cout << "initial precision bits = " << initial_precision_bits << ", iterations required = " << iterations_required << std::endl;
        // initial precision bits = 9, iterations required = 0
        } // if
        else
        {
          iterations_required = 4; // Worst cases should not need more than 4 Halley iterations.
        }

        RealType w0;
        if ( (x > (std::numeric_limits<float>::max)())
           // If x > FLOAT_MAX, then too big to fit into a float, so can't use lambert_w_approx<float>.
          || (x > 1 / epsilon)
          // Multiprecision type where FLOAT_MAX ~=1e38 is too small to hold 1/epsilon.
          )
        { // then can't use lambert_w_approx<float>, so use alternative approximation log(x) - log(log(x).
          w0 = lambert_w_ln_approx(x);
        }
        else
        {
#ifndef BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS
          // RealType argument type can be converted to float.
          w0 = lambert_w_approx(static_cast<float>(x));
#else
          // Not convertible so must use lambert_w_approx<RealType>(x)
          w0 = lambert_w_approx(x);
#endif
        } //

        if (!boost::math::isfinite(w0))
        { // 1st approximation is not finite (unexpected), so quit.
          return policies::raise_domain_error<RealType>(
            "boost::math::lambert_w<%1%>",
            "The parameter %1% approximation is NaN.",
            x,
            Policy());
          return std::numeric_limits<RealType>::quiet_NaN();
        }
        if (iterations_required == 0)
        { // Just return approximation without any refinement.
          return w0;
        }
      RealType w1 = 0; // Refined estimate.
      // TODO only = 0 to avoid  warning C4701: potentially uninitialized local variable 'w1' used.
      RealType previous_diff = boost::math::tools::max_value<RealType>(); // Record the previous diff so can confirm improvement.
      int iterations = 0;
      RealType tolerance = 1 * std::numeric_limits<RealType>::epsilon(); // Or computed get_epsilon??

      while ((iterations < iterations_required) && (iterations <= 10)) // Absolute limit 10 during testing - looping if need this.
      { // Iterate a few times to refine value using Halley's method.
        // Inline Halley iteration, rather than calling boost::math::tools::halley_iterate
        // since we can simplify the expressions algebraically,
        // and don't need most of the error checking of the boilerplate version
        // as we know in advance that the function is reasonably well-behaved,
        // and in any case the derivatives require evaluation of Lambert W!

        RealType expw0 = exp(w0); // Compute from best Lambert W estimate so far.
        // Hope that  w0 * expw0 == x;
        RealType diff = (w0 * expw0) - x; // Difference from x.
        iterations++;
        if (iterations > max_iterations)
        { // Record the maximum iterations used by all calls of Lambert_w - only used during performance testing.
          max_iterations = iterations;
        }
        total_iterations++; // Record the total iterations used by all calls of Lambert_w - only used during performance testing.
        // Halley's method from Luu equation 6.39, line 17.
        // https://en.wikipedia.org/wiki/Halley%27s_method
        // http://www.wolframalpha.com/input/?i=differentiate+productlog(x)
        // differentiate productlog(x)  d/dx(W(x)) = (W(x))/(x W(x) + x) = e^w (1 + w)
        // Wolfram Alpha (d^2 )/(dw^2)(w exp(w) - z) = e^w (w + 2)
        // f'(w) = e^w (1 + w)
        // f''(w) = e^w (2 + w) ,
        // f''(w)/f'(w) = (2+w) / (1+w),  Luu equation 6.38.

        w1 = w0 // Refine a new estimate using Halley's method using Luu equation 6.39.
          - diff / ((expw0 * (w0 + 1) - (w0 + 2) * diff / (w0 + w0 + 2)));
      #ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
        std::cout.precision(std::numeric_limits<RealType>::max_digits10);
        std::cout << "Iteration #" << iterations << ", w0 " << w0 << ", w1 = " << w1;

        if (diff != static_cast<RealType>(0))
        {
          std::cout << ", difference = " << w1 - w0
            << ", relative " << ((w0 / w1) - static_cast<RealType>(1))
            << ", float distance = " << abs(float_distance((w0 * expw0), x)) << std::endl;
        }
        else
        {
          std::cout << ", exact." << std::endl;
        }
        std::cout << "f'(x) = " << diff / (expw0 * (w0 + 1)) << ", f''(x) = " << -diff / ((expw0 * (w0 + 1) - (w0 + 2) * diff / (w0 + w0 + 2))) << std::endl;
      #endif // BOOST_MATH_INSTRUMENT_LAMBERT_W

        if ( // Reached estimate of Lambert W within relative tolerance (usually an epsilon or few).
          (fabs((w0 / w1) - static_cast<RealType>(1)) < tolerance)
          || // Or latest estimate is not better, or worse, so avoid oscillating result (common when near singularity).
          (abs(diff) >= abs(previous_diff))
          )
        {
      #ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
          std::cout << "Return refined within tolerance " << w1 << " after " << iterations << " iterations";
          if (diff != static_cast<RealType>(0))
          {
            std::cout << ", difference = " << ((w0 / w1) - static_cast<RealType>(1));
            std::cout << ", Float distance = " << boost::math::float_distance<RealType>(w0, w1) << std::endl;
          }
          else
          {
            std::cout << ", exact." << std::endl;
          }
      #endif

          if (diff > worst_diff)
          {
            worst_diff = static_cast<double>(diff);
          }
          total_diffs += static_cast<double>(diff);

          return w1; // As good as it can be.
        } // if reached tolerance.

        w0 = w1;
        previous_diff = diff; // Remember so that can check if any new estimate is better.
        }  // while

#ifdef BOOST_MATH_INSTRUMENT_LAMBERT_W
        std::cout << "Not within tolerance " << w1 << " after " << iterations << " iterations" << ", difference = " << previous_diff << std::endl;
#endif
        return w1;
      } // lambert_w_imp
    } // namespace detail

      // User functions.

    // User-defined policy.
    template <class T, class Policy>
    inline typename tools::promote_args<T>::type
    lambert_w(T z, const Policy& pol)
    {
      return detail::lambert_w_imp(z, pol);
    }

    // Default policy.
    template <class T>
    inline typename
    tools::promote_args<T>::type  lambert_w(T z)
    {
      return lambert_w(z, policies::policy<>());
    }

  } // namespace math
} // namespace boost

#endif // BOOST_MATH_SF_LAMBERT_W_HPP
