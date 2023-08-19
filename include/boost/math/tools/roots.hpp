//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_NEWTON_SOLVER_HPP
#define BOOST_MATH_TOOLS_NEWTON_SOLVER_HPP

#ifdef _MSC_VER
#pragma once
#endif
#include <boost/math/tools/complex.hpp> // test for multiprecision types in complex Newton

#include <utility>
#include <cmath>
#include <tuple>
#include <cstdint>

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/cxx03_warn.hpp>
#include <boost/math/tools/throw_exception.hpp>

#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/policies/error_handling.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
using boost::multiprecision::float128;
#endif

namespace boost {
namespace math {
namespace tools {

namespace detail {

namespace dummy {

   template<int n, class T>
   typename T::value_type get(const T&) BOOST_MATH_NOEXCEPT(T);
}

template <class Tuple, class T>
void unpack_tuple(const Tuple& t, T& a, T& b) BOOST_MATH_NOEXCEPT(T)
{
   using dummy::get;
   // Use ADL to find the right overload for get:
   a = get<0>(t);
   b = get<1>(t);
}
template <class Tuple, class T>
void unpack_tuple(const Tuple& t, T& a, T& b, T& c) BOOST_MATH_NOEXCEPT(T)
{
   using dummy::get;
   // Use ADL to find the right overload for get:
   a = get<0>(t);
   b = get<1>(t);
   c = get<2>(t);
}

template <class Tuple, class T>
inline void unpack_0(const Tuple& t, T& val) BOOST_MATH_NOEXCEPT(T)
{
   using dummy::get;
   // Rely on ADL to find the correct overload of get:
   val = get<0>(t);
}

template <class T, class U, class V>
inline void unpack_tuple(const std::pair<T, U>& p, V& a, V& b) BOOST_MATH_NOEXCEPT(T)
{
   a = p.first;
   b = p.second;
}
template <class T, class U, class V>
inline void unpack_0(const std::pair<T, U>& p, V& a) BOOST_MATH_NOEXCEPT(T)
{
   a = p.first;
}

template <class F, class T>
void handle_zero_derivative(F f,
   T& last_f0,
   const T& f0,
   T& delta,
   T& result,
   T& guess,
   const T& min,
   const T& max) noexcept(BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   if (last_f0 == 0)
   {
      // this must be the first iteration, pretend that we had a
      // previous one at either min or max:
      if (result == min)
      {
         guess = max;
      }
      else
      {
         guess = min;
      }
      unpack_0(f(guess), last_f0);
      delta = guess - result;
   }
   if (sign(last_f0) * sign(f0) < 0)
   {
      // we've crossed over so move in opposite direction to last step:
      if (delta < 0)
      {
         delta = (result - min) / 2;
      }
      else
      {
         delta = (result - max) / 2;
      }
   }
   else
   {
      // move in same direction as last step:
      if (delta < 0)
      {
         delta = (result - max) / 2;
      }
      else
      {
         delta = (result - min) / 2;
      }
   }
}

} // namespace

template <class F, class T, class Tol, class Policy>
std::pair<T, T> bisect(F f, T min, T max, Tol tol, std::uintmax_t& max_iter, const Policy& pol) noexcept(policies::is_noexcept_error_policy<Policy>::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   T fmin = f(min);
   T fmax = f(max);
   if (fmin == 0)
   {
      max_iter = 2;
      return std::make_pair(min, min);
   }
   if (fmax == 0)
   {
      max_iter = 2;
      return std::make_pair(max, max);
   }

   //
   // Error checking:
   //
   static const char* function = "boost::math::tools::bisect<%1%>";
   if (min >= max)
   {
      return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function,
         "Arguments in wrong order in boost::math::tools::bisect (first arg=%1%)", min, pol));
   }
   if (fmin * fmax >= 0)
   {
      return boost::math::detail::pair_from_single(policies::raise_evaluation_error(function,
         "No change of sign in boost::math::tools::bisect, either there is no root to find, or there are multiple roots in the interval (f(min) = %1%).", fmin, pol));
   }

   //
   // Three function invocations so far:
   //
   std::uintmax_t count = max_iter;
   if (count < 3)
      count = 0;
   else
      count -= 3;

   while (count && (0 == tol(min, max)))
   {
      T mid = (min + max) / 2;
      T fmid = f(mid);
      if ((mid == max) || (mid == min))
         break;
      if (fmid == 0)
      {
         min = max = mid;
         break;
      }
      else if (sign(fmid) * sign(fmin) < 0)
      {
         max = mid;
      }
      else
      {
         min = mid;
         fmin = fmid;
      }
      --count;
   }

   max_iter -= count;

#ifdef BOOST_MATH_INSTRUMENT
   std::cout << "Bisection required " << max_iter << " iterations.\n";
#endif

   return std::make_pair(min, max);
}

template <class F, class T, class Tol>
inline std::pair<T, T> bisect(F f, T min, T max, Tol tol, std::uintmax_t& max_iter)  noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   return bisect(f, min, max, tol, max_iter, policies::policy<>());
}

template <class F, class T, class Tol>
inline std::pair<T, T> bisect(F f, T min, T max, Tol tol) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::uintmax_t m = (std::numeric_limits<std::uintmax_t>::max)();
   return bisect(f, min, max, tol, m, policies::policy<>());
}


////// Motivation for the Bisection namespace. //////
//
// What's the best way to bisect between a lower bound (lb) and an upper
// bound (ub) during root finding? Let's consider options...
//
// Arithmetic bisection:
//   - The natural choice, but it doesn't always work well. For example, if
//     lb = 1.0 and ub = std::numeric_limits<float>::max(), many bisections
//     may be needed to converge if the root is near 1.
//
// Geometric bisection:
//   - This approach performs much better for the example above, but it
//     too has issues. For example, if lb = 0.0, it gets stuck at 0.0.
//     It also fails if lb and ub have different signs.
//
// In addition to the limitations outlined above, neither of these approaches
// works if ub is infinity. We want a more robust way to handle bisection
// for general root finding problems. That's what this namespace is for.
//
namespace detail {
namespace Bisection {

   ////// The Midpoint754 class //////
   //
   // On a conceptual level, this class is designed to solve the following root
   // finding problem.
   //   - A function f(x) has a single root x_solution somewhere in the interval
   //     [-infinity, +infinity]. For all values below x_solution f(x) is -1. 
   //     For all values above x_solution f(x) is +1. The best way to root find
   //     this problem is to bisect in bit space.
   //
   // Efficient bit space bisection is possible because of the IEEE 754 standard.
   // According to the standard, the bits in floating point numbers are partitioned
   // into three parts: sign, exponent, and mantissa. As long as the sign of the
   // of the number stays the same, increasing numbers in bit space have increasing
   // floating point values starting at zero, and ending at infinity! The table
   // below shows select numbers for float (single precision).
   //
   // 0            |  0 00000000 00000000000000000000000  |  positive zero
   // 1.4013e-45   |  0 00000000 00000000000000000000001  |  std::numeric_limits<float>::denorm_min()
   // 1.17549e-38  |  0 00000001 00000000000000000000000  |  std::numeric_limits<float>::min()
   // 1.19209e-07  |  0 01101000 00000000000000000000000  |  std::numeric_limits<float>::epsilon()
   // 1            |  0 01111111 00000000000000000000000  |  positive one
   // 3.40282e+38  |  0 11111110 11111111111111111111111  |  std::numeric_limits<float>::max()
   // inf          |  0 11111111 00000000000000000000000  |  std::numeric_limits<float>::infinity()
   // nan          |  0 11111111 10000000000000000000000  |  std::numeric_limits<float>::quiet_NaN()
   // 
   // Negative values are similar, but the sign bit is set to 1. By keeping track of the possible
   // sign flip, it can bisect numbers with different signs.
   //
   template <typename T, typename U>
   class Midpoint754 {
   private:
      // Does the bisection in bit space for IEEE 754 floating point numbers.
      // Infinities are allowed. It's assumed that neither x nor y is NaN.
      static_assert(std::numeric_limits<T>::is_iec559, "Type must be IEEE 754 floating point.");
      static_assert(std::is_unsigned<U>::value, "U must be an unsigned integer type.");
      static_assert(sizeof(T) == sizeof(U), "Type and uint size must be the same.");

      // The sign bit is 1, all other bits are 0
      static constexpr U sign_mask() { return U(1) << (sizeof(U) * 8 - 1); }

      // Convert float to bits
      static U float_to_uint(T x) {
         U bits;
         std::memcpy(&bits, &x, sizeof(U));
         return bits;
      }

      // Convert bits to float
      static T uint_to_float(U bits) {
         T x;
         std::memcpy(&x, &bits, sizeof(T));
         return x;
      }

      // Avoids linking against libquadmath for __float128
      static T fabs_imp(T x) {
         U bits_x = float_to_uint(x);
         bits_x = bits_x & ~sign_mask();  // Zero the sign bit
         return uint_to_float(bits_x);
      }

   public:
      static T solve(T x, T y) {
         // Sort so that y has the larger magnitude
         if (fabs_imp(y) < fabs_imp(x)) {
            std::swap(x, y);
         }

         // Recast as unsigned integers
         const U bits_x = float_to_uint(x);
         const U bits_y = float_to_uint(y);

         // Get just the sign bit
         const U sign_x = bits_x & sign_mask();
         const U sign_y = bits_y & sign_mask();

         // Get everything but the sign bit
         const U mag_x = bits_x & ~sign_mask();
         const U mag_y = bits_y & ~sign_mask();

         // Calculate the average magnitude in bits
         U bits_mag = (sign_x == sign_y) ? (mag_y + mag_x) : (mag_y - mag_x);
         bits_mag = bits_mag >> 1;  // Divide by 2

         bits_mag = bits_mag | sign_y;  // Make the sign of bits_mag match the sign of y

         return uint_to_float(bits_mag);
      }

      // NOTE: boost::multiprecision::float128 is cast to __float128
      template <typename V>
      static V solve(V x, V y) { return solve(static_cast<T>(x), static_cast<T>(y)); }

      // Must evaluate to true in order to bisect correctly with infinity
      // Ideally this should be a static assert.
      static bool is_one_plus_max_bits_inf() {
         const U bits_max = float_to_uint((std::numeric_limits<T>::max)());
         const U bits_one_plus_max = bits_max + 1;
         return uint_to_float(bits_one_plus_max) == std::numeric_limits<T>::infinity();
      }

      using type_float = T;  // Used for unit tests
   };  // class Midpoint754


   template <typename T>
   class MidpointNon754 {
   private:
      // NOTE: The Midpoint754 solver should be used when possible because it is faster
      //       than this solver. The two criteria below must be satifsied to use the Midpoint754
      //       solver:
      //           1. The type T must conform to the IEEE 754 standard and
      //           2. The following sequence of steps must produce `numeric_limits<T>::infinity`.
      //              Start with `numeric_limits<T>::max()`. Then reinterpret this value as an
      //              unsigned integer. Next add one to this unsigned integer. Finally reinterpret
      //              the result as type T. This result must equal infinity.
      //       The above two criteria are true for the following datatypes: `float`, `double`,
      //       and `__float128`. The 80 bit variation of `long double` does not satisfy the
      //       second criteria.
      static_assert(!std::is_same<T, float>::value, "Need to use Midpoint754 solver when T is float");
      static_assert(!std::is_same<T, double>::value, "Need to use Midpoint754 solver when T is double");
#if defined(BOOST_HAS_INT128) && defined(BOOST_HAS_FLOAT128)
      static_assert(!std::is_same<T, __float128>::value, "Need to use Midpoint754 solver when T is __float128");
#endif

   public:
      static T solve(T x, T y) {
         const T sign_x = sign(x);
         const T sign_y = sign(y);
         
         // Sign flip return zero
         if (sign_x * sign_y == -1) { return T(0.0); }

         // At least one is positive
         if (0 < sign_x + sign_y) { return do_solve(x, y); }

         // At least one is negative
         return -do_solve(-x, -y);
      }

   private:      
      struct EqZero {
         EqZero(T x) { BOOST_MATH_ASSERT_MSG(x == 0, "x must be zero.");}
      };

      struct EqInf {
         // NOTE: (x == infinity) is checked this way to support types (e.g., boost::math::concepts::real_concept)
         //       that have infinity, but for which this the query std::numeric_limits<T>::infinity() returns
         //       0 due to the std::numeric_limits interface not being implemented.
         EqInf(T x) { BOOST_MATH_ASSERT_MSG(0 < x && !(boost::math::isfinite)(x), "x must be positive infinity."); }
      };

      class PosFinite {
      public:
         PosFinite(T x) : x_(x) {
            BOOST_MATH_ASSERT_MSG(0 < x, "x must be positive.");
            BOOST_MATH_ASSERT_MSG((boost::math::isfinite)(x), "x must be finite.");
         }

         T value() const { return x_; }

      private:
         T x_;
      };

      // Two unknowns
      static T do_solve(T x, T y) {
         if (y < x) {
            return do_solve(y, x);
         }

         if (x == 0) {
            return do_solve(EqZero(x), y);  // Zero and ???
         } else if ((boost::math::isfinite)(x)) {
            return do_solve(PosFinite(x), y);  // Finite and ???
         } else {
            return x;  // Infinity and infinity
         }
      }

      // One unknowns
      static T do_solve(EqZero x, T y) {
         if (y == 0) {
            return T(0.0);  // Zero and zero
         } else if ((boost::math::isfinite)(y)) {
            return do_solve(x, PosFinite(y));  // Zero and finite
         } else {
            return T(1.5);  // Zero and infinity
         }
      }
      static T do_solve(PosFinite x, T y) {
         if ((boost::math::isfinite)(y)) {
            return do_solve(x, PosFinite(y));  // Finite and finite
         } else {
            return do_solve(x, EqInf(y));  // Finite and infinity
         }
      }

      // Zero unknowns
      template <typename U = T>
      static typename std::enable_if<std::numeric_limits<U>::is_specialized, T>::type
      do_solve(PosFinite x, EqInf y) {
          return do_solve(x, PosFinite((std::numeric_limits<U>::max)()));
      }
      template <typename U = T>
      static typename std::enable_if<!std::numeric_limits<U>::is_specialized, T>::type
      do_solve(PosFinite x, EqInf y) {
         BOOST_MATH_THROW_EXCEPTION(std::runtime_error("infinite bounds support requires specialization"));
         return T(0); // Unreachable, but suppresses warnings
      }

      template <typename U = T>
      static typename std::enable_if<std::numeric_limits<U>::is_specialized, U>::type
      do_solve(EqZero x, PosFinite y) { 
         const auto get_smallest_value = []() {
            const U denorm_min = std::numeric_limits<U>::denorm_min();
            if (denorm_min != 0) { return denorm_min; }

            // NOTE: the two lines below can be removed when
            //       https://github.com/boostorg/multiprecision/pull/562
            //       gets merged.
            const U min = (std::numeric_limits<U>::min)();
            if (min != 0) { return min; }

            BOOST_MATH_THROW_EXCEPTION(std::runtime_error("This type probably isn't actually specialized."));
            return T(0); // Unreachable, but suppresses warnings
         };

         return do_solve(PosFinite(get_smallest_value()), y);
      }
      template <typename U = T>
      static typename std::enable_if<!std::numeric_limits<U>::is_specialized, U>::type
      do_solve(EqZero x, PosFinite y) {
         // This function quickly gets a value that is small relative to y.
         const auto fn_appx_denorm_min = [](U y){
            T accumulator = T(0.5);
            for (int i = 1; i < 100; ++i) {
               if (y * accumulator == 0) {
                  return y;
               } else {
                  y = y * accumulator;
                  accumulator = accumulator * T(0.5);
               }
            }
            return y;
         };

         const U denorm_min = fn_appx_denorm_min(y.value());
         return do_solve(PosFinite(denorm_min), y);
      }

      static T do_solve(PosFinite x, PosFinite y) {
         BOOST_MATH_ASSERT_MSG(x.value() <= y.value(), "x must be less than or equal to y.");

         const T value_x = x.value();
         const T value_y = y.value();
         
         // Take arithmetic mean if they are close enough
         if (value_y < value_x * 8) { return (value_y - value_x) / 2 + value_x; }  // NOTE: avoids overflow

         // Take geometric mean if they are far apart
         using std::sqrt;
         return sqrt(value_x) * sqrt(value_y);  // NOTE: avoids overflow
      }
   }; // class MidpointNon754

   // The purposes of this class is to not cause compiler warnings from unused functions.
   class CalcMidpoint {
   public:
      template <typename T>
      static MidpointNon754<T> get_solver(T) {
         return MidpointNon754<T>();
      }
      static Midpoint754<float, std::uint32_t> get_solver(float) {
         return Midpoint754<float, std::uint32_t>();
      }
      static Midpoint754<double, std::uint64_t> get_solver(double) {
         return Midpoint754<double, std::uint64_t>();
      }
#if defined(BOOST_HAS_INT128) && defined(BOOST_HAS_FLOAT128)
      static Midpoint754<__float128, boost::uint128_type> get_solver(__float128) {
         return Midpoint754<__float128, boost::uint128_type>();
      }
      static Midpoint754<__float128, boost::uint128_type> get_solver(boost::multiprecision::float128) {
         return Midpoint754<__float128, boost::uint128_type>();
      }
#endif
      template <typename T>
      static T calc_midpoint(T x, T y) {
         return get_solver(x).solve(x, y);
      }
   };

}  // namespace Bisection
}  // namespace detail

template <class F, class T>
T newton_raphson_iterate(F f, T guess, T min, T max, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   BOOST_MATH_STD_USING

   static const char* function = "boost::math::tools::newton_raphson_iterate<%1%>";
   if (min > max)
   {
      return policies::raise_evaluation_error(function, "Range arguments in wrong order in boost::math::tools::newton_raphson_iterate(first arg=%1%)", min, boost::math::policies::policy<>());
   }

   T f0(0), f1, last_f0(0);
   T result = guess;

   T factor = static_cast<T>(ldexp(1.0, 1 - digits));
   T delta = tools::max_value<T>();
   T delta1 = tools::max_value<T>();
   T delta2 = tools::max_value<T>();

   //
   // We use these to sanity check that we do actually bracket a root,
   // we update these to the function value when we update the endpoints
   // of the range.  Then, provided at some point we update both endpoints
   // checking that max_range_f * min_range_f <= 0 verifies there is a root
   // to be found somewhere.  Note that if there is no root, and we approach 
   // a local minima, then the derivative will go to zero, and hence the next
   // step will jump out of bounds (or at least past the minima), so this
   // check *should* happen in pathological cases.
   //
   T max_range_f = 0;
   T min_range_f = 0;

   std::uintmax_t count(max_iter);

#ifdef BOOST_MATH_INSTRUMENT
   std::cout << "Newton_raphson_iterate, guess = " << guess << ", min = " << min << ", max = " << max
      << ", digits = " << digits << ", max_iter = " << max_iter << "\n";
#endif

   do {
      last_f0 = f0;
      delta2 = delta1;
      delta1 = delta;
      if (count == 0) {
         return policies::raise_evaluation_error(function, "Ran out of iterations in boost::math::tools::newton_raphson_iterate, guess: %1%", guess, boost::math::policies::policy<>());
      } else {
         --count;
      }
      detail::unpack_tuple(f(result), f0, f1);

      if (0 == f0)
         break;
      if (f1 == 0)
      {
         // Oops zero derivative!!!
         detail::handle_zero_derivative(f, last_f0, f0, delta, result, guess, min, max);
      }
      else
      {
         delta = f0 / f1;
      }
#ifdef BOOST_MATH_INSTRUMENT
      std::cout << "Newton iteration " << max_iter - count << ", delta = " << delta << ", residual = " << f0 << "\n";
#endif
      if (fabs(delta * 2) > fabs(delta2))
      {
         // Last two steps haven't converged.
         const T x_other = (delta > 0) ? min : max;
         delta = result - detail::Bisection::CalcMidpoint::calc_midpoint(result, x_other);
         // reset delta1/2 so we don't take this branch next time round:
         delta1 = 3 * delta;
         delta2 = 3 * delta;
      }
      guess = result;
      result -= delta;
      if (result <= min)
      {
         delta = 0.5F * (guess - min);
         result = guess - delta;
         if ((result == min) || (result == max))
            break;
      }
      else if (result >= max)
      {
         delta = 0.5F * (guess - max);
         result = guess - delta;
         if ((result == min) || (result == max))
            break;
      }
      // Update brackets:
      if (delta > 0)
      {
         max = guess;
         max_range_f = f0;
      }
      else if (delta < 0)  // Cannot have "else" here, as delta being zero is not indicative of failure
      {
         min = guess;
         min_range_f = f0;
      }
      //
      // Sanity check that we bracket the root:
      //
      if (max_range_f * min_range_f > 0)
      {
         return policies::raise_evaluation_error(function, "There appears to be no root to be found in boost::math::tools::newton_raphson_iterate, perhaps we have a local minima near current best guess of %1%", guess, boost::math::policies::policy<>());
      }
   } while (fabs(result * factor) < fabs(delta) || result == 0);

   max_iter -= count;

#ifdef BOOST_MATH_INSTRUMENT
   std::cout << "Newton Raphson required " << max_iter << " iterations\n";
#endif

   return result;
}

template <class F, class T>
inline T newton_raphson_iterate(F f, T guess, T min, T max, int digits) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::uintmax_t m = (std::numeric_limits<std::uintmax_t>::max)();
   return newton_raphson_iterate(f, guess, min, max, digits, m);
}

namespace detail {

   struct halley_step
   {
      template <class T>
      static T step(const T& /*x*/, const T& f0, const T& f1, const T& f2) noexcept(BOOST_MATH_IS_FLOAT(T))
      {
         using std::fabs;
         T denom = 2 * f0;
         T num = 2 * f1 - f0 * (f2 / f1);
         T delta;

         BOOST_MATH_INSTRUMENT_VARIABLE(denom);
         BOOST_MATH_INSTRUMENT_VARIABLE(num);

         if ((fabs(num) < 1) && (fabs(denom) >= fabs(num) * tools::max_value<T>()))
         {
            // possible overflow, use Newton step:
            delta = f0 / f1;
         }
         else
            delta = denom / num;
         return delta;
      }
   };

   template <class F, class T>
   T bracket_root_towards_min(F f, T guess, const T& f0, T& min, T& max, std::uintmax_t& count) noexcept(BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())));

   template <class F, class T>
   T bracket_root_towards_max(F f, T guess, const T& f0, T& min, T& max, std::uintmax_t& count) noexcept(BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
   {
      using std::fabs;
      using std::ldexp;
      using std::abs;
      using std::frexp;
      if(count < 2)
         return guess - (max + min) / 2; // Not enough counts left to do anything!!
      //
      // Move guess towards max until we bracket the root, updating min and max as we go:
      //
      int e;
      frexp(max / guess, &e);
      e = abs(e);
      T guess0 = guess;
      T multiplier = e < 64 ? static_cast<T>(2) : static_cast<T>(ldexp(T(1), e / 32));
      T f_current = f0;
      if (fabs(min) < fabs(max))
      {
         while (--count && ((f_current < 0) == (f0 < 0)))
         {
            min = guess;
            guess *= multiplier;
            if (guess > max)
            {
               guess = max;
               f_current = -f_current;  // There must be a change of sign!
               break;
            }
            multiplier *= e > 1024 ? 8 : 2;
            unpack_0(f(guess), f_current);
         }
      }
      else
      {
         //
         // If min and max are negative we have to divide to head towards max:
         //
         while (--count && ((f_current < 0) == (f0 < 0)))
         {
            min = guess;
            guess /= multiplier;
            if (guess > max)
            {
               guess = max;
               f_current = -f_current;  // There must be a change of sign!
               break;
            }
            multiplier *= e > 1024 ? 8 : 2;
            unpack_0(f(guess), f_current);
         }
      }

      if (count)
      {
         max = guess;
         if (multiplier > 16)
            return (guess0 - guess) + bracket_root_towards_min(f, guess, f_current, min, max, count);
      }
      return guess0 - (max + min) / 2;
   }

   template <class F, class T>
   T bracket_root_towards_min(F f, T guess, const T& f0, T& min, T& max, std::uintmax_t& count) noexcept(BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
   {
      using std::fabs;
      using std::ldexp;
      using std::abs;
      using std::frexp;
      if (count < 2)
         return guess - (max + min) / 2; // Not enough counts left to do anything!!
      //
      // Move guess towards min until we bracket the root, updating min and max as we go:
      //
      int e;
      frexp(guess / min, &e);
      e = abs(e);
      T guess0 = guess;
      T multiplier = e < 64 ? static_cast<T>(2) : static_cast<T>(ldexp(T(1), e / 32));
      T f_current = f0;

      if (fabs(min) < fabs(max))
      {
         while (--count && ((f_current < 0) == (f0 < 0)))
         {
            max = guess;
            guess /= multiplier;
            if (guess < min)
            {
               guess = min;
               f_current = -f_current;  // There must be a change of sign!
               break;
            }
            multiplier *= e > 1024 ? 8 : 2;
            unpack_0(f(guess), f_current);
         }
      }
      else
      {
         //
         // If min and max are negative we have to multiply to head towards min:
         //
         while (--count && ((f_current < 0) == (f0 < 0)))
         {
            max = guess;
            guess *= multiplier;
            if (guess < min)
            {
               guess = min;
               f_current = -f_current;  // There must be a change of sign!
               break;
            }
            multiplier *= e > 1024 ? 8 : 2;
            unpack_0(f(guess), f_current);
         }
      }

      if (count)
      {
         min = guess;
         if (multiplier > 16)
            return (guess0 - guess) + bracket_root_towards_max(f, guess, f_current, min, max, count);
      }
      return guess0 - (max + min) / 2;
   }

   template <class Stepper, class F, class T>
   T second_order_root_finder(F f, T guess, T min, T max, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
   {
      BOOST_MATH_STD_USING

#ifdef BOOST_MATH_INSTRUMENT
        std::cout << "Second order root iteration, guess = " << guess << ", min = " << min << ", max = " << max
        << ", digits = " << digits << ", max_iter = " << max_iter << "\n";
#endif
      static const char* function = "boost::math::tools::halley_iterate<%1%>";
      if (min >= max)
      {
         return policies::raise_evaluation_error(function, "Range arguments in wrong order in boost::math::tools::halley_iterate(first arg=%1%)", min, boost::math::policies::policy<>());
      }

      T f0(0), f1, f2;
      T result = guess;

      T factor = ldexp(static_cast<T>(1.0), 1 - digits);
      T delta = (std::max)(T(10000000 * guess), T(10000000));  // arbitrarily large delta
      T last_f0 = 0;
      T delta1 = delta;
      T delta2 = delta;
      bool out_of_bounds_sentry = false;

   #ifdef BOOST_MATH_INSTRUMENT
      std::cout << "Second order root iteration, limit = " << factor << "\n";
   #endif

      //
      // We use these to sanity check that we do actually bracket a root,
      // we update these to the function value when we update the endpoints
      // of the range.  Then, provided at some point we update both endpoints
      // checking that max_range_f * min_range_f <= 0 verifies there is a root
      // to be found somewhere.  Note that if there is no root, and we approach 
      // a local minima, then the derivative will go to zero, and hence the next
      // step will jump out of bounds (or at least past the minima), so this
      // check *should* happen in pathological cases.
      //
      T max_range_f = 0;
      T min_range_f = 0;

      std::uintmax_t count(max_iter);

      do {
         last_f0 = f0;
         delta2 = delta1;
         delta1 = delta;
#ifndef BOOST_NO_EXCEPTIONS
         try
#endif
         {
            detail::unpack_tuple(f(result), f0, f1, f2);
         }
#ifndef BOOST_NO_EXCEPTIONS
         catch (const std::overflow_error&)
         {
            f0 = max > 0 ? tools::max_value<T>() : -tools::min_value<T>();
            f1 = f2 = 0;
         }
#endif
         --count;

         BOOST_MATH_INSTRUMENT_VARIABLE(f0);
         BOOST_MATH_INSTRUMENT_VARIABLE(f1);
         BOOST_MATH_INSTRUMENT_VARIABLE(f2);

         if (0 == f0)
            break;
         if (f1 == 0)
         {
            // Oops zero derivative!!!
            detail::handle_zero_derivative(f, last_f0, f0, delta, result, guess, min, max);
         }
         else
         {
            if (f2 != 0)
            {
               delta = Stepper::step(result, f0, f1, f2);
               if (delta * f1 / f0 < 0)
               {
                  // Oh dear, we have a problem as Newton and Halley steps
                  // disagree about which way we should move.  Probably
                  // there is cancelation error in the calculation of the
                  // Halley step, or else the derivatives are so small
                  // that their values are basically trash.  We will move
                  // in the direction indicated by a Newton step, but
                  // by no more than twice the current guess value, otherwise
                  // we can jump way out of bounds if we're not careful.
                  // See https://svn.boost.org/trac/boost/ticket/8314.
                  delta = f0 / f1;
                  if (fabs(delta) > 2 * fabs(result))
                     delta = (delta < 0 ? -1 : 1) * 2 * fabs(result);
               }
            }
            else
               delta = f0 / f1;
         }
   #ifdef BOOST_MATH_INSTRUMENT
         std::cout << "Second order root iteration, delta = " << delta << ", residual = " << f0 << "\n";
   #endif
         T convergence = fabs(delta / delta2);
         if ((convergence > 0.8) && (convergence < 2))
         {
            // last two steps haven't converged.
            if (fabs(min) < 1 ? fabs(1000 * min) < fabs(max) : fabs(max / min) > 1000)
            {
               if(delta > 0)
                  delta = bracket_root_towards_min(f, result, f0, min, max, count);
               else
                  delta = bracket_root_towards_max(f, result, f0, min, max, count);
            }
            else
            {
               delta = (delta > 0) ? (result - min) / 2 : (result - max) / 2;
               if ((result != 0) && (fabs(delta) > result))
                  delta = sign(delta) * fabs(result) * 0.9f; // protect against huge jumps!
            }
            // reset delta2 so that this branch will *not* be taken on the
            // next iteration:
            delta2 = delta * 3;
            delta1 = delta * 3;
            BOOST_MATH_INSTRUMENT_VARIABLE(delta);
         }
         guess = result;
         result -= delta;
         BOOST_MATH_INSTRUMENT_VARIABLE(result);

         // check for out of bounds step:
         if (result < min)
         {
            T diff = ((fabs(min) < 1) && (fabs(result) > 1) && (tools::max_value<T>() / fabs(result) < fabs(min)))
               ? T(1000)
               : (fabs(min) < 1) && (fabs(tools::max_value<T>() * min) < fabs(result))
               ? ((min < 0) != (result < 0)) ? -tools::max_value<T>() : tools::max_value<T>() : T(result / min);
            if (fabs(diff) < 1)
               diff = 1 / diff;
            if (!out_of_bounds_sentry && (diff > 0) && (diff < 3))
            {
               // Only a small out of bounds step, lets assume that the result
               // is probably approximately at min:
               delta = 0.99f * (guess - min);
               result = guess - delta;
               out_of_bounds_sentry = true; // only take this branch once!
            }
            else
            {
               if (fabs(float_distance(min, max)) < 2)
               {
                  result = guess = (min + max) / 2;
                  break;
               }
               delta = bracket_root_towards_min(f, guess, f0, min, max, count);
               result = guess - delta;
               if (result <= min)
                  result = float_next(min);
               if (result >= max)
                  result = float_prior(max);
               guess = min;
               continue;
            }
         }
         else if (result > max)
         {
            T diff = ((fabs(max) < 1) && (fabs(result) > 1) && (tools::max_value<T>() / fabs(result) < fabs(max))) ? T(1000) : T(result / max);
            if (fabs(diff) < 1)
               diff = 1 / diff;
            if (!out_of_bounds_sentry && (diff > 0) && (diff < 3))
            {
               // Only a small out of bounds step, lets assume that the result
               // is probably approximately at min:
               delta = 0.99f * (guess - max);
               result = guess - delta;
               out_of_bounds_sentry = true; // only take this branch once!
            }
            else
            {
               if (fabs(float_distance(min, max)) < 2)
               {
                  result = guess = (min + max) / 2;
                  break;
               }
               delta = bracket_root_towards_max(f, guess, f0, min, max, count);
               result = guess - delta;
               if (result >= max)
                  result = float_prior(max);
               if (result <= min)
                  result = float_next(min);
               guess = min;
               continue;
            }
         }
         // update brackets:
         if (delta > 0)
         {
            max = guess;
            max_range_f = f0;
         }
         else
         {
            min = guess;
            min_range_f = f0;
         }
         //
         // Sanity check that we bracket the root:
         //
         if (max_range_f * min_range_f > 0)
         {
            return policies::raise_evaluation_error(function, "There appears to be no root to be found in boost::math::tools::newton_raphson_iterate, perhaps we have a local minima near current best guess of %1%", guess, boost::math::policies::policy<>());
         }
      } while(count && (fabs(result * factor) < fabs(delta)));

      max_iter -= count;

   #ifdef BOOST_MATH_INSTRUMENT
      std::cout << "Second order root finder required " << max_iter << " iterations.\n";
   #endif

      return result;
   }
} // T second_order_root_finder

template <class F, class T>
T halley_iterate(F f, T guess, T min, T max, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   return detail::second_order_root_finder<detail::halley_step>(f, guess, min, max, digits, max_iter);
}

template <class F, class T>
inline T halley_iterate(F f, T guess, T min, T max, int digits) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::uintmax_t m = (std::numeric_limits<std::uintmax_t>::max)();
   return halley_iterate(f, guess, min, max, digits, m);
}

namespace detail {

   struct schroder_stepper
   {
      template <class T>
      static T step(const T& x, const T& f0, const T& f1, const T& f2) noexcept(BOOST_MATH_IS_FLOAT(T))
      {
         using std::fabs;
         T ratio = f0 / f1;
         T delta;
         if ((x != 0) && (fabs(ratio / x) < 0.1))
         {
            delta = ratio + (f2 / (2 * f1)) * ratio * ratio;
            // check second derivative doesn't over compensate:
            if (delta * ratio < 0)
               delta = ratio;
         }
         else
            delta = ratio;  // fall back to Newton iteration.
         return delta;
      }
   };

}

template <class F, class T>
T schroder_iterate(F f, T guess, T min, T max, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   return detail::second_order_root_finder<detail::schroder_stepper>(f, guess, min, max, digits, max_iter);
}

template <class F, class T>
inline T schroder_iterate(F f, T guess, T min, T max, int digits) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::uintmax_t m = (std::numeric_limits<std::uintmax_t>::max)();
   return schroder_iterate(f, guess, min, max, digits, m);
}
//
// These two are the old spelling of this function, retained for backwards compatibility just in case:
//
template <class F, class T>
T schroeder_iterate(F f, T guess, T min, T max, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   return detail::second_order_root_finder<detail::schroder_stepper>(f, guess, min, max, digits, max_iter);
}

template <class F, class T>
inline T schroeder_iterate(F f, T guess, T min, T max, int digits) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::uintmax_t m = (std::numeric_limits<std::uintmax_t>::max)();
   return schroder_iterate(f, guess, min, max, digits, m);
}

#ifndef BOOST_NO_CXX11_AUTO_DECLARATIONS
/*
   * Why do we set the default maximum number of iterations to the number of digits in the type?
   * Because for double roots, the number of digits increases linearly with the number of iterations,
   * so this default should recover full precision even in this somewhat pathological case.
   * For isolated roots, the problem is so rapidly convergent that this doesn't matter at all.
   */
template<class Complex, class F>
Complex complex_newton(F g, Complex guess, int max_iterations = std::numeric_limits<typename Complex::value_type>::digits)
{
   typedef typename Complex::value_type Real;
   using std::norm;
   using std::abs;
   using std::max;
   // z0, z1, and z2 cannot be the same, in case we immediately need to resort to Muller's Method:
   Complex z0 = guess + Complex(1, 0);
   Complex z1 = guess + Complex(0, 1);
   Complex z2 = guess;

   do {
      auto pair = g(z2);
      if (norm(pair.second) == 0)
      {
         // Muller's method. Notation follows Numerical Recipes, 9.5.2:
         Complex q = (z2 - z1) / (z1 - z0);
         auto P0 = g(z0);
         auto P1 = g(z1);
         Complex qp1 = static_cast<Complex>(1) + q;
         Complex A = q * (pair.first - qp1 * P1.first + q * P0.first);

         Complex B = (static_cast<Complex>(2) * q + static_cast<Complex>(1)) * pair.first - qp1 * qp1 * P1.first + q * q * P0.first;
         Complex C = qp1 * pair.first;
         Complex rad = sqrt(B * B - static_cast<Complex>(4) * A * C);
         Complex denom1 = B + rad;
         Complex denom2 = B - rad;
         Complex correction = (z1 - z2) * static_cast<Complex>(2) * C;
         if (norm(denom1) > norm(denom2))
         {
            correction /= denom1;
         }
         else
         {
            correction /= denom2;
         }

         z0 = z1;
         z1 = z2;
         z2 = z2 + correction;
      }
      else
      {
         z0 = z1;
         z1 = z2;
         z2 = z2 - (pair.first / pair.second);
      }

      // See: https://math.stackexchange.com/questions/3017766/constructing-newton-iteration-converging-to-non-root
      // If f' is continuous, then convergence of x_n -> x* implies f(x*) = 0.
      // This condition approximates this convergence condition by requiring three consecutive iterates to be clustered.
      Real tol = (max)(abs(z2) * std::numeric_limits<Real>::epsilon(), std::numeric_limits<Real>::epsilon());
      bool real_close = abs(z0.real() - z1.real()) < tol && abs(z0.real() - z2.real()) < tol && abs(z1.real() - z2.real()) < tol;
      bool imag_close = abs(z0.imag() - z1.imag()) < tol && abs(z0.imag() - z2.imag()) < tol && abs(z1.imag() - z2.imag()) < tol;
      if (real_close && imag_close)
      {
         return z2;
      }

   } while (max_iterations--);

   // The idea is that if we can get abs(f) < eps, we should, but if we go through all these iterations
   // and abs(f) < sqrt(eps), then roundoff error simply does not allow that we can evaluate f to < eps
   // This is somewhat awkward as it isn't scale invariant, but using the Daubechies coefficient example code,
   // I found this condition generates correct roots, whereas the scale invariant condition discussed here:
   // https://scicomp.stackexchange.com/questions/30597/defining-a-condition-number-and-termination-criteria-for-newtons-method
   // allows nonroots to be passed off as roots.
   auto pair = g(z2);
   if (abs(pair.first) < sqrt(std::numeric_limits<Real>::epsilon()))
   {
      return z2;
   }

   return { std::numeric_limits<Real>::quiet_NaN(),
            std::numeric_limits<Real>::quiet_NaN() };
}
#endif


#if !defined(BOOST_NO_CXX17_IF_CONSTEXPR)
// https://stackoverflow.com/questions/48979861/numerically-stable-method-for-solving-quadratic-equations/50065711
namespace detail
{
#if defined(BOOST_GNU_STDLIB) && !defined(_GLIBCXX_USE_C99_MATH_TR1)
inline float fma_workaround(float x, float y, float z) { return ::fmaf(x, y, z); }
inline double fma_workaround(double x, double y, double z) { return ::fma(x, y, z); }
#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline long double fma_workaround(long double x, long double y, long double z) { return ::fmal(x, y, z); }
#endif
#endif            
template<class T>
inline T discriminant(T const& a, T const& b, T const& c)
{
   T w = 4 * a * c;
#if defined(BOOST_GNU_STDLIB) && !defined(_GLIBCXX_USE_C99_MATH_TR1)
   T e = fma_workaround(-c, 4 * a, w);
   T f = fma_workaround(b, b, -w);
#else
   T e = std::fma(-c, 4 * a, w);
   T f = std::fma(b, b, -w);
#endif
   return f + e;
}

template<class T>
std::pair<T, T> quadratic_roots_imp(T const& a, T const& b, T const& c)
{
#if defined(BOOST_GNU_STDLIB) && !defined(_GLIBCXX_USE_C99_MATH_TR1)
   using boost::math::copysign;
#else
   using std::copysign;
#endif
   using std::sqrt;
   if constexpr (std::is_floating_point<T>::value)
   {
      T nan = std::numeric_limits<T>::quiet_NaN();
      if (a == 0)
      {
         if (b == 0 && c != 0)
         {
            return std::pair<T, T>(nan, nan);
         }
         else if (b == 0 && c == 0)
         {
            return std::pair<T, T>(0, 0);
         }
         return std::pair<T, T>(-c / b, -c / b);
      }
      if (b == 0)
      {
         T x0_sq = -c / a;
         if (x0_sq < 0) {
            return std::pair<T, T>(nan, nan);
         }
         T x0 = sqrt(x0_sq);
         return std::pair<T, T>(-x0, x0);
      }
      T discriminant = detail::discriminant(a, b, c);
      // Is there a sane way to flush very small negative values to zero?
      // If there is I don't know of it.
      if (discriminant < 0)
      {
         return std::pair<T, T>(nan, nan);
      }
      T q = -(b + copysign(sqrt(discriminant), b)) / T(2);
      T x0 = q / a;
      T x1 = c / q;
      if (x0 < x1)
      {
         return std::pair<T, T>(x0, x1);
      }
      return std::pair<T, T>(x1, x0);
   }
   else if constexpr (boost::math::tools::is_complex_type<T>::value)
   {
      typename T::value_type nan = std::numeric_limits<typename T::value_type>::quiet_NaN();
      if (a.real() == 0 && a.imag() == 0)
      {
         using std::norm;
         if (b.real() == 0 && b.imag() && norm(c) != 0)
         {
            return std::pair<T, T>({ nan, nan }, { nan, nan });
         }
         else if (b.real() == 0 && b.imag() && c.real() == 0 && c.imag() == 0)
         {
            return std::pair<T, T>({ 0,0 }, { 0,0 });
         }
         return std::pair<T, T>(-c / b, -c / b);
      }
      if (b.real() == 0 && b.imag() == 0)
      {
         T x0_sq = -c / a;
         T x0 = sqrt(x0_sq);
         return std::pair<T, T>(-x0, x0);
      }
      // There's no fma for complex types:
      T discriminant = b * b - T(4) * a * c;
      T q = -(b + sqrt(discriminant)) / T(2);
      return std::pair<T, T>(q / a, c / q);
   }
   else // Most likely the type is a boost.multiprecision.
   {    //There is no fma for multiprecision, and in addition it doesn't seem to be useful, so revert to the naive computation.
      T nan = std::numeric_limits<T>::quiet_NaN();
      if (a == 0)
      {
         if (b == 0 && c != 0)
         {
            return std::pair<T, T>(nan, nan);
         }
         else if (b == 0 && c == 0)
         {
            return std::pair<T, T>(0, 0);
         }
         return std::pair<T, T>(-c / b, -c / b);
      }
      if (b == 0)
      {
         T x0_sq = -c / a;
         if (x0_sq < 0) {
            return std::pair<T, T>(nan, nan);
         }
         T x0 = sqrt(x0_sq);
         return std::pair<T, T>(-x0, x0);
      }
      T discriminant = b * b - 4 * a * c;
      if (discriminant < 0)
      {
         return std::pair<T, T>(nan, nan);
      }
      T q = -(b + copysign(sqrt(discriminant), b)) / T(2);
      T x0 = q / a;
      T x1 = c / q;
      if (x0 < x1)
      {
         return std::pair<T, T>(x0, x1);
      }
      return std::pair<T, T>(x1, x0);
   }
}
}  // namespace detail

template<class T1, class T2 = T1, class T3 = T1>
inline std::pair<typename tools::promote_args<T1, T2, T3>::type, typename tools::promote_args<T1, T2, T3>::type> quadratic_roots(T1 const& a, T2 const& b, T3 const& c)
{
   typedef typename tools::promote_args<T1, T2, T3>::type value_type;
   return detail::quadratic_roots_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(c));
}

#endif

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_NEWTON_SOLVER_HPP
