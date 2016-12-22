// Copyright Thomas Luu, 2014
// Copyright Paul A. Bristow 2016

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

//#define BOOST_MATH_INSTRUMENT  // #define for diagnostic output.

#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/constants/constants.hpp>

#include <boost/fixed_point/fixed_point.hpp> // fixed_point

#define NEWTON

#ifdef BOOST_MATH_INSTRUMENT
#  include <iostream>  // Only for testing.
#endif
#include <limits> // numeric_limits Inf Nan ...
#include <cmath> // exp
#include <type_traits> // is_integral


// Define a temporary constant for exp(-1) for use in checks at the singularity.
namespace boost
{ namespace math 
  {
    namespace constants
    { // constant exp(-1) = 0.367879441 is where x = -exp(-1) is a singularity.
      BOOST_DEFINE_MATH_CONSTANT(expminusone, 3.67879441171442321595523770161460867e-01, "0.36787944117144232159552377016146086744581113103176783450783680169746149574489980335714727434591964374662732527");
    }
  } 
}


namespace boost {
  namespace fixed_point {

    //! \brief Specialization for log1p for fixed_point.

    /*! \details This simple (and fast) approximation gives a result within a few epsilon/ULP of the nearest representable value.
     This makes it suitable for use as a first approximation by the Lambert W function below.
     It only uses (and thus relies on for accuracy of) the standard log function.
     Refinement of the Lambert W estimate using Halley's method
     seems to get quickly to high accuracy in a very few iterations,
     so that small inaccuracy in the 1st appproximation are unimportant.
     It avoid any Taylor series refinements that involve high powers of x
     that would easily go outside the often limited range of fixed_point types.
    */

    template<const int IntegralRange, const int FractionalResolution, typename RoundMode, typename OverflowMode>
    negatable<IntegralRange, FractionalResolution, RoundMode, OverflowMode>
      log1p(negatable<IntegralRange, FractionalResolution, RoundMode, OverflowMode> x)
    {
      // https://github.com/f32c/arduino/blob/master/hardware/fpga/f32c/system/src/math/log1p.c
      // HP - 15C Advanced Functions Handbook, p.193.

      typedef negatable<IntegralRange, FractionalResolution, RoundMode, OverflowMode> local_negatable_type;

      local_negatable_type unity(1);
      // A fixed_point type night not include unity.  Does something sensible happen?
      local_negatable_type u = unity + x;
      if (u == unity)
      { //
        return x;
      }
      else
      {
        local_negatable_type result = log(u) * (x / (u - unity));
#ifdef BOOST_MATH_INSTRUMENT
        std::cout << "log1p fp approx " << result << std::endl;
#endif
        // For example: x = 0.5, HP 0.78843689, diff -1.52587891e-05
        return result;
      } // if else
    } // log1p_fp

  } //namespace fixed_point
} // namespace boost


namespace local
{
  // log1p_dodgy version seems to work OK, but...
  /*  J M cautions that "The formula relies on:

  1) That the operations are not re-ordered.

  2) That each operation is correctly rounded.
    Is true only for floating-point arithmetic unless you compile with
    special flags to force full-IEEE compatible mode,
    but that slows the whole program down...
   Not a concern for fixed-point?

  3) That compiler optimisations don't perform symbolic math on the expressions and simplify them."
  */

template<typename T>
T log1p_dodgy(T x)
{
  // log(1+x) == (log(1+x) * x) / ((1-x) - 1)
  T result = -(log(1 + x) * x) / ((1 - x) - 1); // From note in Boost.Math log1p

  // \boost\libs\math\include\boost\math\special_functions\log1p.hpp
  // Algorithm log1p is part of C99, but is not yet provided by many compilers.
  //
  // The Boost.Math version uses a Taylor series expansion for 0.5 > x > epsilon, which may
  // require up to std::numeric_limits<T>::digits+1 terms to be calculated.
  // It would be much more efficient to use the equivalence:
  //   log(1+x) == (log(1+x) * x) / ((1-x) - 1)
  // Unfortunately many optimizing compilers make such a mess of this, that
  // it performs no better than log(1+x): which is to say not very well at all.

  return result;
} // template<typename T> T log1p(T x)

template<typename T>
T log1p(T x)
{
  //! \brief log1p function.
  //! This is probably faster than using boost::math::log1p (if a little less accurate)
  //! and should be good enough for computing initial estimate.
  //! \details Formula from a Note in
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
  //

  // It might be possible to use Newton or Halley's method to refine this approximation?
  // But improvement may not be useful for estimating the Lambert W value because it is close enough
  // for the Halley's method to converge just as well as with a better initial estimate.

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
    double dx = static_cast<double>(x);
    double lp1x = log1p(dx);
#ifdef BOOST_MATH_INSTRUMENT
    std::cout << "true " << lp1x
      << ", log1p_fp " << log1p(x)
      //<< ", dodgy " << dodgy << ", HP " << result
      //<< ", diff " << dodgy - result << "epsilons = " << (dodgy - result)/std::numeric_limits<T>::epsilon()
      << std::endl;
    // For example: x = 0.5, dodgy 0.788421631, HP 0.78843689, diff -1.52587891e-05
#endif

    return result;
  }
}  // template<typename T> T log1p(T x)

} // namespace local

namespace boost
{
  namespace math
  {
    //! \brief Lambert W function.
    //! \details Based on Thomas Luu, Thesis (2016).
    //! Notes refer to equations, page 96-98 and C++ code from algorithm Routine 11.

    template <class RealType = double, class Policy = policies::policy<> >
    RealType productlog(RealType x)
    {
      BOOST_MATH_STD_USING  // for ADL of std functions.

      using boost::math::log1p;
      using boost::math::constants::root_two; // 1.414 ...
      using boost::math::constants::e; // e(1) = 2.71828...
      using boost::math::constants::one_div_root_two; // 0.707106
      using boost::math::constants::expminusone; // 0.36787944

      std::cout.precision(std::numeric_limits <RealType>::max_digits10);
      std::cout << std::showpoint << std::endl; // Show all trailing zeros.

      // Catch the very common mistake of providing an integer value as parameter to productlog.
      // need to ensure it is a floating-point type (of the desired type, float 1.f, double 1., or long double 1.L),
      // or static_cast, for example:  static_cast<float>(1) or static_cast<cpp_dec_float_50>(1).
      // Want to allow fixed_point types too.
      BOOST_STATIC_ASSERT_MSG(!std::is_integral<RealType>::value, "Must be floating-point, not integer type, for example W(1.), not W(1)!");

#ifdef BOOST_MATH_INSTRUMENT
      std::cout << std::showpoint << std::endl; // Show all trailing zeros.
      //std::cout.precision(std::numeric_limits<RealType>::max_digits10); // Show all possibly significant digits.
      std::cout.precision(std::numeric_limits<RealType>::digits10); // Show all significant digits.
#endif
      // Check on range of x.

      //std::cout << "-exp(-1) = " << -expminusone<RealType>() << std::endl;

      // Special case of -exp(-1)) // -0.3678794411714423215955237701614608674458111310
      // Can't use if (x < -exp(-1)) because 1-bit difference in accuracy of exp means is inconsistent.
      // if (x < static_cast<RealType>(-0.3678794411714423215955237701614608674458111310L) )
      if (x < -expminusone<RealType>())
      { // If x <  -0.367879 then W(x) would be complex (not handled with this implementation).
        // Or might throw an exception?  Use domain error for policy here.

        //std::cout << "Would be Complex " << x << ' ' << static_cast<RealType>(-0.3678794411714423215955237701614608674458111310L) << " return NaN!" << std::endl;
        return std::numeric_limits<RealType>::quiet_NaN();
      }
      //else if (x == static_cast<RealType>(-0.3678794411714423215955237701614608674458111310))
      else if (x == -expminusone<RealType>() )
      {
        //std::cout << "At Singularity " << x << ' ' << static_cast<RealType>(-0.3678794411714423215955237701614608674458111310L) << " returned " << static_cast<RealType>(-1) << std::endl;

        return static_cast<RealType>(-1);
      }

      if (x == 0)
      { // Special case of zero (for speed and to avoid log(0)and /0 etc).
        // TODO check that values very close to zero do not fail?
        return static_cast<RealType>(0);
      }

      RealType w0;
      // Upper branch W0 is divided in two regions: -1 <= w0_minus <=0 and 0 <= w0_plus.
      if (x > 0)
      { // Luu Routine 11, line 8, and equation 6.44, from Barry et al.
        // (1.2 and 2.4 are 'exact' from integer fractions 6/5 and 12/5).
        w0 = log(static_cast<RealType>(1.2L) * x / log(static_cast<RealType>(2.4L) * x / log1p(static_cast<RealType>(2.4L) * x)));
        //w0 = log(static_cast<RealType>(1.2L) * x / log(static_cast<RealType>(2.4L) * x / local::log1p(static_cast<RealType>(2.4L) * x)));
        // local::log1p does seem to refine OK with multiprecision.
      }
      else
      { // > -exp(-1) or > -0.367879
        // so for real result need x > -0.367879 >= 0
      // Might treat near -exp(-1) == -0.367879 values differently as approximations may not quite converge?

        RealType sqrt_v = root_two<RealType>() * sqrt(static_cast<RealType>(1) + e<RealType>() * x); // nu = sqrt(2 + 2e * z) Line 10.

        RealType n2 = static_cast<RealType>(3) * root_two<RealType>() + static_cast<RealType>(6); // 3 * sqrt(2) + 6, Line 11.
        // == 10.242640687119285146

        RealType n2num = (static_cast<RealType>(2237) + (static_cast<RealType>(1457) * root_two<RealType>())) * e<RealType>()
          - (static_cast<RealType>(4108) * root_two<RealType>()) - static_cast<RealType>(5764); //  Numerator Line 11.

        RealType n2denom = (static_cast<RealType>(215) + (static_cast<RealType>(199) * root_two<RealType>())) * e<RealType>()
          - (430 * root_two<RealType>()) - static_cast<RealType>(796); // Denominator Line 11.
        RealType nn2 = n2num / n2denom; // Line 11 full fraction.
        nn2 *= sqrt_v; // Line 11 last part.
        n2 -= nn2; //  Line 11 complete, equation 6.44, from Barry et al (erratum).
        RealType n1 = (static_cast<RealType>(1) - one_div_root_two<RealType>()) * (n2 + root_two<RealType>()); // Line 12.

        w0 = -1 + sqrt_v * (n2 + sqrt_v) / (n2 + sqrt_v + n1 * sqrt_v); // W0 1st approximation, Line 13, equation 6.40.
      }

      #ifdef BOOST_MATH_INSTRUMENT
        // std::cout << "\n1st approximation  = " << w0 << std::endl;
      #endif

      int iterations = 0;
      RealType tolerance = 3 * std::numeric_limits<RealType>::epsilon();
      RealType w1;

      while (iterations <= 5)
      { // Iterate a few times to refine value using Halley's method.
        RealType expw0 = exp(w0); // Compute from best W estimate so far.
        // Hope that  w0 * expw0 == x;
        RealType diff = w0 * expw0 - x; // Difference from x.
        if (abs(diff) <= tolerance/4)
        { // Too close for Halley iteration to improve?
          break; // Avoid oscillating forever around value.
        }

        // Newton/Raphson's method https://en.wikipedia.org/wiki/Newton%27s_method
        // f(w) = w e^w -z = 0  Luu equation 6.37
        // f'(w) = e^w (1 + w), Wolfram alpha (d)/(dw)(f(w) = w exp(w) - z) = e^w (w + 1)
        // f(w) / f'(w)
        // w1 = w0 - (expw0 * (w0 + 1)); // Refine new Newton/Raphson estimate.
        // Takes typically 6 iterations to converge within tolerance,
        // whereas Halley usually takes 1 to 3 iterations,
        // so Newton unlikely to be quicker than additional computation cost of 2nd derivative.
        // Might use Newton if near to x ~= -exp(-1) == -0.367879


        // Halley's method from Luu equation 6.39, line 17.
        // https://en.wikipedia.org/wiki/Halley%27s_method
        // f''(w) = e^w (2 + w) , Wolfram Alpha (d^2 )/(dw^2)(w exp(w) - z) = e^w (w + 2)
        // f''(w) / f'(w) = (2+w) / (1+w),  Luu equation 6.38.

        w1 = w0 // Refine new Halley estimate.
            - diff /
            ((expw0 * (w0 + 1) - (w0 + 2) * diff / (w0 + w0 + 2))); // Luu equation 6.39.

        if (fabs((w0 / w1) - 1) < tolerance)
        { // Reached estimate of Lambert W within tolerance (usually an epsilon or few).
          break;
        }

#ifdef BOOST_MATH_INSTRUMENT
        std::cout.precision(std::numeric_limits<RealType>::digits10);
        std::cout <<"Iteration #" << iterations << ", w0 " << w0 << ", w1 = " << w1 << ", difference = " << diff << ", relative " << (w0 / w1 - static_cast<RealType>(1)) << std::endl;
        std::cout << "f'(x) = " << diff / (expw0 * (w0 + 1)) << ", f''(x) = " << - diff / ((expw0 * (w0 + 1) - (w0 + 2) * diff / (w0 + w0 + 2))) << std::endl;
        //std::cout << "Newton new = " << wn << std::endl;
#endif
        w0 = w1;
        iterations++;
      } // while

#ifdef BOOST_MATH_INSTRUMENT
      std::cout << "Final " << w1 << " after " << iterations << " iterations" << ", difference = " << (w0 / w1 - static_cast<RealType>(1)) << std::endl;
#endif
      return w1;
    }

  } // namespace math
} // namespace boost
