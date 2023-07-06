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

#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/policies/error_handling.hpp>

#include <boost/math/tools/minima.hpp>
#include <iostream>

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


#if 0

   //
   template <class C>
   class Optional {
   public:
      Optional() : has_data_(false) {}
      Optional(C c) : has_data_(true), c_(c) {}

      bool has_data() const { return has_data_; }
      const C& c() const { 
         assert(has_data_);
         return c_;
      }

      void replace(C c) { c_ = c; }

   private:
      bool has_data_;
      C c_;
   };

   // Enum for bracketing states: case_1, case_2U, case_2B
   enum class BracketState { case_1, case_2U, case_2B, case_2I };

   // Enum for solve outputs
   enum class OutputResult { success, failure };

   // Solve output type
   template <class T>
   using SolveOutput = std::pair<OutputResult, T>;






#endif













namespace detail {


#if 1

   // constexpr bool DEBUG_PRINT = false;
   constexpr bool DEBUG_PRINT = true;


   // Stores the last two Newton step sizes
   template <class T>
   class StepSizeHistory {
      public:
         StepSizeHistory() { reset(); }

         void reset() {
            // NOTE: for T = boost::math::concepts::real_concept, both std::numeric_limits<T>::infinity()
            // and std::numeric_limits<T>::max() are 0
            dx_p_ = static_cast<T>(1) / 0;
            dx_pp_ = static_cast<T>(1) / 0;
         }

         void push(T input) {
            dx_pp_ = dx_p_;
            dx_p_ = input;
         }

         T dx_p() const { return dx_p_; }
         T dx_pp() const { return dx_pp_; }

      private:
         T dx_p_;
         T dx_pp_;
   };

   // Stores x, f(x), and f'(x)
   template <class T, class Step>
   class CacheOrder1 {
   public:
      template <class F>
      CacheOrder1(F f, T x) : x_(x), has_f0_f1_(true) { detail::unpack_tuple(f(x_), f0_, f1_); }
      CacheOrder1(T x) : x_(x) , has_f0_f1_(false) {}

      T calc_dx() const { return Step().calc_dx(*this); }

      bool is_l() const { return 0 <= calc_dx(); }
      bool is_h() const { return !is_l(); }
 
      bool is_B(const CacheOrder1<T, Step>& z) const { return !is_same_sign(z); }
      bool is_I(const CacheOrder1<T, Step>& z) const { return !is_B(z) && !is_U(z); }
      bool is_U(const CacheOrder1<T, Step>& z) const {
         if (is_ordered(z)) {
            return is_U_ordered(*this, z);
         } else {
            return is_U_ordered(z, *this);
         }
      }
      std::pair<CacheOrder1<T, Step>, CacheOrder1<T, Step>> sort(const CacheOrder1<T, Step>& z) const {
         if (is_ordered(z)) {
            return std::make_pair(*this, z);
         } else {
            return std::make_pair(z, *this);
         }
      }
      bool is_ordered(const CacheOrder1<T, Step>& z) const { return x() < z.x(); }
      bool is_U_ordered(const CacheOrder1<T, Step>& l, const CacheOrder1<T, Step>& h) const {
         const bool is_pos = 0 < l.f0();

         T l_slope = l.f1();
         T h_slope = h.f1();

         if (!is_pos) {
            l_slope *= -1;
            h_slope *= -1;
         }

         return (l_slope < 0) && (0 < h_slope);
      }
      bool is_smaller(const CacheOrder1<T, Step>& z) const {
         using std::fabs;
         return fabs(f0_) < fabs(z.f0());
      }

      bool is_same_sign(const CacheOrder1<T, Step>& z) const { return sign(f0_) == sign(z.f0()); }
      bool is_same_slope(const CacheOrder1<T, Step>& z) const { return sign(f1_) == sign(z.f1()); }

      void print() const {
         if (!DEBUG_PRINT) return;

         std::cout << "x: " << x_;
         if (has_f0_f1_) {
            assert(has_f0_f1_);
            std::cout << ", f0: " << f0_ << ", f1: " << f1_;
         }
         std::cout << std::endl;
      }

      bool is_x_only() const { return !has_f0_f1_; }
      bool is_x_f0_f1() const { return has_f0_f1_; }

      bool is_dx_small(const T& factor) const {
         // |dx| ≤ |x| * factor
         // |f0| / |f1| ≤ |x| * factor
         // |f0| ≤ |x| * |f1| * factor
         using std::fabs;

         const bool bool_eval   = fabs(f0_) <= fabs(x_ * f1_) * factor;
         const bool bool_eval_2 = fabs(f0_) <= fabs(x_ * f1_ * factor);

         if (DEBUG_PRINT) {
            std::cout.precision(200);
            std::cout << "fabs(f0_):                   " << fabs(f0_) << std::endl;
            std::cout << "fabs(x_ * f1_):              " << fabs(x_ * f1_) << std::endl;
            std::cout << "factor:                      " << factor << std::endl;
            std::cout << "bool_eval;                   " << bool_eval << std::endl;
            std::cout << "bool_eval_2;                 " << bool_eval_2 << std::endl;
            std::cout << "fabs(f0_) == 0               " << (fabs(f0_) == 0) << std::endl;
            std::cout << "fabs(x_ * f1_) * factor == 0 " << (fabs(x_ * f1_) * factor == 0) << std::endl;

            // Print data type T to std::cout
            std::cout << "T: " << typeid(T).name() << std::endl;
            std::cout << std::endl;
         }


         return bool_eval;
      }

      T x() const { return x_; }
      T f0() const { assert(has_f0_f1_); return f0_; }  // f(x)
      T f1() const { assert(has_f0_f1_); return f1_; }  // f'(x)

   private:
      T x_;
      T f0_;  // f(x)
      T f1_;  // f'(x)
      bool has_f0_f1_;
   };

   //
   //
   template <class T, class Step>
   class NEXT {
   public:
      NEXT<T, Step> static FromX(T x) {
         return NEXT<T, Step>(x);
      }

      NEXT<T, Step> static TakeStep(const CacheOrder1<T, Step>& cx) {
         const T dx = cx.calc_dx();
         return NEXT<T, Step>(cx.x() + dx, dx);
      }

      NEXT<T, Step> static From_x_dx(T x, T dx) {
         return NEXT<T, Step>(x, dx);
      }

      T x() const { return x_; }
      T dx() const { return dx_; }

      void print() const {
         if (!DEBUG_PRINT) return;
         std::cout << "X_DX_NEXT -- x: " << x_;
         std::cout << ", dx: " << dx_;
         std::cout << std::endl;
      }

   private:
      NEXT(T x) : x_(x), dx_(static_cast<T>(1) / 0) {}
      NEXT(T x, T dx) : x_(x), dx_(dx) {}
      
      T x_;
      T dx_;
   };

   // Calculates the step with Newton's method
   template <class T>
   class StepNewton {
   public:
      T calc_dx(const CacheOrder1<T, StepNewton>& z) const { return -z.f0() / z.f1(); }
   };

   // Enum for status flags
   enum class StatusFlag { normal, dx_small, out_of_iterations, 
      case_1, case_2U, case_2B, case_2I,
      success, failure, stall, f0_is_zero, case_0
   };



////////////////////////////////////////////////////////////////////////////////////////


   //
   //
   template <class F, class T, class Step>
   class ContextF {
   public:
      ContextF(F f, T xl, T xh, int digits, std::uintmax_t max_iter)
         : f_(f)
         , factor_(static_cast<T>(ldexp(1.0, 1 - digits)))
         , max_iter_(max_iter)
         , count_(max_iter)
         , l_orig_(xl)
         , h_orig_(xh)
      {
         max_iter_ = std::min(max_iter_, static_cast<std::uintmax_t>(400));
         count_ = std::min(count_, static_cast<std::uintmax_t>(400));
      }

      CacheOrder1<T, Step> operator()(const NEXT<T, Step>& next) {
         const T x = next.x();

         if (count_ <= 0) {
            static const char* function = "boost::math::tools::detail::ContextF<%1%>";
            policies::raise_evaluation_error(function, "Unable to bracket root in boost::math::tools::detail::ContextF, last closest guess was %1%", x, boost::math::policies::policy<>());
         }

         if (DEBUG_PRINT) std::cout << "EVAL IS: " << count_ << ", at: " << x <<  std::endl;

         --count_;
         return CacheOrder1<T, Step>(f_, x);
      }

      T dx_p() const { return dx_history_.dx_p(); }
      T dx_pp() const { return dx_history_.dx_pp(); }
      std::uintmax_t iters_used() const { return max_iter_ - count_; }

      T factor() const { return factor_; }
      std::uintmax_t max_iter() const { return max_iter_; }
      std::uintmax_t count() const { return count_; }
      StepSizeHistory<T>& dx_history() { return dx_history_; }

   private:
      F f_;
      T factor_;
      std::uintmax_t max_iter_;
      std::uintmax_t count_;
      StepSizeHistory<T> dx_history_;
      T l_orig_;
      T h_orig_;
   };



////////////////////////////////////////////////////////////////////////////////////////

   //
   //
   template <class F, class T, class Step>
   class BoundState {
   public:
      BoundState(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) 
         : cache_l_(xl)
         , cache_h_(xh)
         , flag_(StatusFlag::case_1)
         , context_(f, xl, xh, digits, max_iter)
      {
         const auto next = NEXT<T, Step>::FromX(x);
         const auto cache_x_temp = context_(next);
         update_cache_sign_dx(cache_x_temp);
      }

      T calc_dx_true(const CacheOrder1<T, Step>& cx) const {
         return cx.x() - cache_x().x();
      }

      // CacheOrder1<T, Step> operator()(T x) { return context_(x); }
      CacheOrder1<T, Step> operator()(const NEXT<T, Step>& next) { return context_(next); }
      
      void print_dx_history() const {
         if (!DEBUG_PRINT) return;
         // std::cout << "dx_p: " << dx_p() << std::endl;
         // std::cout << "dx_pp: " << dx_pp() << std::endl;
      }

      void update_cache_l(const CacheOrder1<T, Step>& cx) {
         if (DEBUG_PRINT) std::cout << "set dx_true to: " << calc_dx_true(cx) << std::endl;
         if (DEBUG_PRINT) std::cout << "dx_newton:      " << cx.calc_dx() << std::endl;
         dx_history().push(calc_dx_true(cx));
         print_dx_history();

         if (DEBUG_PRINT) std::cout << "UPDATE_CACHE_L to: " << cx.x() << std::endl;

         cache_l_ = cx;
         is_last_eval_l_ = true;
      }
      void update_cache_h(const CacheOrder1<T, Step>& cx) {
         if (DEBUG_PRINT) std::cout << "set dx_true to: " << calc_dx_true(cx) << std::endl;
         if (DEBUG_PRINT) std::cout << "dx_newton:      " << cx.calc_dx() << std::endl;
         dx_history().push(calc_dx_true(cx));
         print_dx_history();

         if (DEBUG_PRINT) std::cout << "UPDATE_CACHE_H to: " << cx.x() << std::endl;

         cache_h_ = cx;
         is_last_eval_l_ = false;
      }

      void update_cache_sign_f0(const CacheOrder1<T, Step>& cx) {
         if (cx.is_same_sign(cache_l_)) {
            update_cache_l(cx);
         } else {
            update_cache_h(cx);
         }
      }

      void update_cache_sign_dx(const CacheOrder1<T, Step>& cx) {
         if (cx.is_l()) {
            update_cache_l(cx);
         } else {
            update_cache_h(cx);
         }
      }

      bool is_converge_too_slow(const T& dx) const {
         using std::fabs;
         return fabs(dx_pp()) < fabs(dx * 4);
      }
      
      T midpoint() const { return (cache_l_.x() + cache_h_.x()) / 2; }
      bool is_inbounds(T x) const { return (l_smart() <= x) && (x <= h_smart()); }

      T l_smart() const { return cache_l_.x(); }
      T h_smart() const { return cache_h_.x(); }

      bool is_dx_small(const NEXT<T, Step>& next) const { 
         using std::fabs;

         const T side_1 = fabs(next.dx());
         const T side_2 = fabs(next.x()) * factor();

         // if (DEBUG_PRINT) std::cout << "side_1: " << side_1 << std::endl;
         // if (DEBUG_PRINT) std::cout << "side_2: " << side_2 << std::endl;

         return side_1 <= side_2;

         // return fabs(dx_p()) <= fabs(cache_x().x()) * factor();
      }

      void set_flag_if_dx_small(const NEXT<T, Step>& next) {
         if (cache_x().f0() == 0) {
            set_flag(StatusFlag::f0_is_zero);
         } else if (is_dx_small(next)) {
            set_flag(StatusFlag::dx_small);
         }
      }
      void set_flag_if_count_zero() {
         if (count() <= 0) { 
            set_flag(StatusFlag::out_of_iterations);
            if (DEBUG_PRINT) std::cout << "OUT OF ITERATIONS: count: " << count() << std::endl;
         }
      }
      void set_flag_if_stall(const NEXT<T, Step>& next) {
         if (next.dx() == 0) { set_flag(StatusFlag::stall); }
      }

      void materialize_l() {
         if (is_last_eval_l_) { return; }
         cache_l_ = (*this)(l_smart());
         is_last_eval_l_ = true;
      }
      void materialize_h() {
         if (!is_last_eval_l_) { return; }
         cache_h_ = (*this)(h_smart());
         is_last_eval_l_ = false;
      }

      const CacheOrder1<T, Step>& cache_x() const {
         if (is_last_eval_l_) { return cache_l_; } else { return cache_h_; }
      }

      int num_x_only() const { return cache_l_.is_x_only() + cache_h_.is_x_only(); }

      T get_limit() const {
         assert(num_x_only() == 1);
         return is_last_eval_l_ ? cache_h_.x() : cache_l_.x();
      }

      void set_state_2() {
         assert(num_x_only() == 0);

         const auto& cache_l = cache_l_;
         const auto& cache_h = cache_h_;

         if (cache_l.is_B(cache_h)) {
            flag_ = StatusFlag::case_2B;
         } else if (cache_l.is_U(cache_h)) {
            flag_ = StatusFlag::case_2U;
         } else {
            assert(cache_l.is_I(cache_h));
            flag_ = StatusFlag::case_2I;
         } 
      }

      void print() const {
         if (!DEBUG_PRINT) return;
         std::cout << "is_last_eval_l_: " << is_last_eval_l_ << std::endl;
         std::cout << "l:               "; cache_l_.print();
         std::cout << "h:               "; cache_h_.print();
         std::cout << "flag:            ";
         print_root_status_flag();
      }

      void print_root_status_flag() const {
         if (!DEBUG_PRINT) return;
         switch (flag()) {
            case StatusFlag::normal: std::cout << "normal" << std::endl; break;
            case StatusFlag::dx_small: std::cout << "dx_small" << std::endl; break;
            case StatusFlag::out_of_iterations: std::cout << "out_of_iterations" << std::endl; break;
            case StatusFlag::case_1: std::cout << "case_1" << std::endl; break;
            case StatusFlag::case_2U: std::cout << "case_2U" << std::endl; break;
            case StatusFlag::case_2B: std::cout << "case_2B" << std::endl; break;
            case StatusFlag::case_2I: std::cout << "case_2I" << std::endl; break;
            case StatusFlag::success: std::cout << "success" << std::endl; break;
            case StatusFlag::failure: std::cout << "failure" << std::endl; break;
            case StatusFlag::stall: std::cout << "stall" << std::endl; break;
            case StatusFlag::f0_is_zero: std::cout << "f0_is_zero" << std::endl; break;
            case StatusFlag::case_0: std::cout << "case_0" << std::endl; break;
         }
      }
      
      void set_flag(StatusFlag flag) { flag_ = flag; }

      bool is_B() const { return cache_l_.is_B(cache_h_); }
      bool is_U() const { return cache_l_.is_U(cache_h_); }
      bool is_I() const { return cache_l_.is_I(cache_h_); }      

      StepSizeHistory<T>& dx_history() { return context_.dx_history(); }
      std::uintmax_t count() const { return context_.count(); }
      T factor() const { return context_.factor(); }
      std::uintmax_t max_iter() const { return context_.max_iter(); }
      T dx_p() const { return context_.dx_p(); }
      T dx_pp() const { return context_.dx_pp(); }
      std::uintmax_t iters_used() const { return context_.iters_used(); }

      NEXT<T, Step> create_next_for_x(T x) {
         const T dx = x - cache_x().x();
         return NEXT<T, Step>::From_x_dx(x, dx);
      }

      const CacheOrder1<T, Step>& cache_l() const { return cache_l_; }
      const CacheOrder1<T, Step>& cache_h() const { return cache_h_; }
      StatusFlag flag() const { return flag_; }

      void set_x_next_value(T x) { x_next_value_ = x; }
      T x_next_value() const { return x_next_value_; }

   private:
      bool is_last_eval_l_;
      CacheOrder1<T, Step> cache_l_;
      CacheOrder1<T, Step> cache_h_;
      StatusFlag flag_;
      ContextF<F, T, Step> context_;
      T x_next_value_;
   };

////////////////////////////////////////////////////////////////////////////////////////

   template <class F, class T, class Step>
   class Bound2B;
   
   template <class F, class T, class Step>
   class Bound2U;

   template <class F, class T, class Step>
   class Bound1;

////////////////////////////////////////////////////////////////////////////////////////

   //
   //
   template <class F, class T, class Step>
   class BoundFn {
   public:
      void do_debug_printout(BoundState<F, T, Step>& b) {
         b.print();
         if (DEBUG_PRINT) std::cout << std::endl;
      }

      void solve_next(BoundState<F, T, Step>& b) {
         do_debug_printout(b);

         while (is_bound_still_valid(b)) {
            const auto next = calc_next(b);
            if (DEBUG_PRINT) next.print(); 

            b.set_x_next_value(next.x());

            if (is_stop(b, next)) { 
               if (DEBUG_PRINT) std::cout << "GOT STOP SIGNAL" << std::endl;
               break;
            }

            const auto cx = b(next);
            this->update_cache(b, cx);
            do_debug_printout(b);
         }
         if (DEBUG_PRINT) std::cout << "DONE WITH LOOP" << std::endl;
      }

      NEXT<T, Step> calc_next(BoundState<F, T, Step>& b) {
         // Attempt Newton step
         const auto& cx = b.cache_x();
         const auto next = NEXT<T, Step>::TakeStep(cx);

         const bool c1 = !(b.is_inbounds(next.x()));
         const bool c2 = b.is_converge_too_slow(next.dx());

         if (c1) {
            if (DEBUG_PRINT) std::cout << "NOT INBOUNDS x_next: " << next.x() << std::endl;
         }
         if (c2) {
            if (DEBUG_PRINT) std::cout << "CONVERGE TOO SLOW dx_next: " << next.dx() << std::endl;
         }

         // if (!(b.is_inbounds(next.x())) || b.is_converge_too_slow(next.dx())) {
         if (c1 || c2) {
            b.dx_history().reset();
            if (DEBUG_PRINT) std::cout << "TOOK BISECTION wanted: " << next.x() << ", dx: " << next.dx() << std::endl;
            return this->calc_next_bisection(b);
         } else {
            if (DEBUG_PRINT) std::cout << "TOOK NEWTON -------------" << std::endl;
            return next;
         }
      }

      NEXT<T, Step> calc_next_bisection_regular(BoundState<F, T, Step>& b) {
         const auto x_mid = b.midpoint();
         if (x_mid == b.l_smart() || x_mid == b.h_smart()) {
            b.set_flag(StatusFlag::dx_small);
         }
         // const auto next = NEXT<T, Step>::FromX(x_mid);
         const auto next = b.create_next_for_x(x_mid);
         return next;
      }

      virtual bool is_stop(BoundState<F, T, Step>& b, const NEXT<T,Step> next) = 0;
      virtual NEXT<T, Step> calc_next_bisection(BoundState<F, T, Step>& b) = 0;
      virtual void update_cache(BoundState<F, T, Step>& b, const CacheOrder1<T, Step>& cx) = 0;
      virtual void print_class() const = 0;
      // virtual void is_bound_still_valid() const = 0;
      virtual bool is_bound_still_valid(const BoundState<F, T, Step>& b) const = 0;
   };

   //
   //
   template <class F, class T, class Step>
   class Bound1 : public BoundFn<F, T, Step> {
   public:
      bool is_stop(BoundState<F, T, Step>& b, const NEXT<T,Step> next) override {
         // if (b.num_x_only() != 2) {
            b.set_flag_if_dx_small(next);
            b.set_flag_if_count_zero();
            b.set_flag_if_stall(next);
         // }
         return b.flag() != StatusFlag::case_1;
      }

      NEXT<T, Step> calc_next_bisection(BoundState<F, T, Step>& b) override {
         const T x = b.cache_x().x();
         const T x_limit = b.get_limit();
         const T dx_limit = x_limit - x;
         const T dx_8 = b.cache_x().calc_dx() * 1000;

         using std::fabs;
         if (fabs(dx_limit) < fabs(dx_8)) {
            return b.create_next_for_x(x_limit);
            // b.set_x_dx_next_for_x(x_limit);
            // b.materialize_l();
            // b.materialize_h();
            // b.set_state_2();
         } else {
            // update_cache(b, b(x + dx_8));
            // b.set_x_dx_next_for_x(x + dx_8);
            return b.create_next_for_x(x + dx_8);
         }
      }

      void update_cache(BoundState<F, T, Step>& b, const CacheOrder1<T, Step>& cx) override {
         b.update_cache_sign_dx(cx);
         if (b.num_x_only() == 0) {
            b.set_state_2();
         } else {
            if (b.get_limit() == cx.x()) {
               b.set_flag(StatusFlag::case_2I);
            }
         }
      }

      void print_class() const override { 
         if (!DEBUG_PRINT) return;
         std::cout << "Bound1" << std::endl; }

      bool is_bound_still_valid(const BoundState<F, T, Step>& b) const override { return b.flag() == StatusFlag::case_1; }
   };

   //
   //
   template <class F, class T, class Step>
   class Bound2B : public BoundFn<F, T, Step> {
   public:
       bool is_stop(BoundState<F, T, Step>& b, const NEXT<T,Step> next) override {
         b.set_flag_if_dx_small(next);
         b.set_flag_if_count_zero();
         b.set_flag_if_stall(next);
         return b.flag() != StatusFlag::case_2B;
      }

      NEXT<T, Step> calc_next_bisection(BoundState<F, T, Step>& b) override {
         return this->calc_next_bisection_regular(b);
      }

      void update_cache(BoundState<F, T, Step>& b, const CacheOrder1<T, Step>& cx) override {
         b.update_cache_sign_f0(cx);
      }

      void print_class() const override { 
         if (!DEBUG_PRINT) return;
         std::cout << "Bound2B" << std::endl; }

      bool is_bound_still_valid(const BoundState<F, T, Step>& b) const override { return b.flag() == StatusFlag::case_2B; }
   };

   //
   //
   template <class F, class T, class Step>
   class Bound2U : public BoundFn<F, T, Step> {
   public:
       bool is_stop(BoundState<F, T, Step>& b, const NEXT<T,Step> next) override {
         b.set_flag_if_dx_small(next);
         b.set_flag_if_count_zero();
         return b.flag() != StatusFlag::case_2U;
      }

      NEXT<T, Step> calc_next_bisection(BoundState<F, T, Step>& b) override {
         return this->calc_next_bisection_regular(b);
      }

      void update_cache(BoundState<F, T, Step>& b, const CacheOrder1<T, Step>& cx) override {
         if (is_dx_bigger_than_prev(b, cx)) {
            b.set_flag(StatusFlag::case_2I);
         }
         b.update_cache_sign_dx(cx);
      }

      bool is_dx_bigger_than_prev(BoundState<F, T, Step>& b, const CacheOrder1<T, Step>& cx) const {
         using std::fabs;
         const auto& cache_same_side = cx.is_l() ? b.cache_l() : b.cache_h();
         const T dx_cx = cx.calc_dx();
         const T dx_same_side = cache_same_side.calc_dx();
         const bool is_dx_bigger = fabs(dx_same_side) < fabs(dx_cx);
         return is_dx_bigger;
      }

      void print_class() const override {
         if (!DEBUG_PRINT) return;
         std::cout << "Bound2U" << std::endl; }
   
      bool is_bound_still_valid(const BoundState<F, T, Step>& b) const override { return b.flag() == StatusFlag::case_2U; }
   };


////////////////
////////////////
////////////////


   // template <class F, class T, class Step>
   // std::pair<bool, T> solve(BoundState<F, T, Step>& state) {

   //    state.is_finished();
   // }

   template <class F, class T, class Step>
   std::pair<bool, T> solve(BoundState<F, T, Step>& state) {

      if (DEBUG_PRINT) {
         std::cout << std::endl;
         std::cout << "=========================================" << std::endl;
         std::cout << "Count: " << state.count() << std::endl;
         std::cout << "STATUS IS: ";
         state.print_root_status_flag();
         std::cout << "-----------------------------------------" << std::endl;
         std::cout << std::endl;
      }

      if (state.flag() == StatusFlag::case_1) {
         auto solver_case_1 = Bound1<F, T, Step>();
         solver_case_1.solve_next(state);
         return solve<F, T, Step>(state);

      } else if (state.flag() == StatusFlag::case_2B) {
         auto solver_case_2b = Bound2B<F, T, Step>();
         solver_case_2b.solve_next(state);
         return solve<F, T, Step>(state);

      } else if (state.flag() == StatusFlag::case_2U) {
         auto solver_case_2u = Bound2U<F, T, Step>();
         solver_case_2u.solve_next(state);
         return solve<F, T, Step>(state);

      } else if (state.flag() == StatusFlag::case_2I) {
         return std::make_pair(false, state.x_next_value());

      } else if (state.flag() == StatusFlag::stall) {
         return std::make_pair(true, state.x_next_value());

      } else if (state.flag() == StatusFlag::dx_small) {
         return std::make_pair(true, state.x_next_value());

      } else if (state.flag() == StatusFlag::f0_is_zero) {
         return std::make_pair(true, state.x_next_value());

      } else if (state.flag() == StatusFlag::out_of_iterations) {
         if (DEBUG_PRINT) std::cout << "ERROR OUT OF ITERATIONS" << std::endl;
         assert(false);

      } else {

         // std::cout << "NOT HANDLED: " << static_cast<int>(state.flag()) << std::endl;
         if (DEBUG_PRINT) std::cout << "NOT HANDLED: " << std::endl;
         state.print_root_status_flag();
         assert(false);
      }

      assert(false);
      return std::make_pair(false, state.x_next_value());
   }



/////////////////////////

   template <class F, class T, class Step>
   T solve_2(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>()))) {
      // Create BoundState
      BoundState<F, T, Step> state(f, x, xl, xh, digits, max_iter);


   }

/////////////////////////



   template <class F, class T, class Step>
   T solve(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>()))) {
      if (DEBUG_PRINT) std::cout << "================= newton_raphson_iterate ==================== with max_iter: " << max_iter << std::endl;
      
      if (DEBUG_PRINT) {
         std::cout << "l_orig: " << xl << std::endl;
         std::cout << "x_orig: " << x << std::endl;
         std::cout << "h_orig: " << xh << std::endl;
      }



      // Create BoundState
#if 1
      BoundState<F, T, Step> state(f, x, xl, xh, digits, max_iter);

#else
      BoundState<F, T, Step> state(f, xl, xh, digits, max_iter);

      // Calc NEXT
      const auto next = NEXT<T, Step>(x);

      // Calc CacheX
      const auto cache_x_orig = state(next);




      // // Set up problem
      // state.update_cache_sign_dx(cache_x_orig);

      // // auto solver_case_1 = Bound1<F, T, Step>();
      // Bound1<F, T, Step>().solve(state);
#endif      











//////////////////////////////////////////////////////////////////////////


      // Get CacheX
      const auto cache_x_orig = state.cache_x();

      // Solve problem
      auto p = solve(state);


//////////////////////////////////////////////////////////////////////////
      // Define lambdas
      auto fn_form_bracket_l = [&](){
         const auto next = NEXT<T, Step>::FromX(xl);
         const auto cache_l_orig = state(next);
         state.update_cache_l(cache_l_orig);
         state.update_cache_h(cache_x_orig);
      };
      auto fn_form_bracket_h = [&](){
         const auto next = NEXT<T, Step>::FromX(xh);
         const auto cache_h_orig = state(next);
         state.update_cache_l(cache_x_orig);
         state.update_cache_h(cache_h_orig);
      };
      auto fn_return = [&](const BoundState<F, T, Step>& s, const std::pair<bool, T> p) {
         s.print_root_status_flag();
         assert(p.first);
         max_iter = s.iters_used();
         return p.second;
      };
//////////////////////////////////////////////////////////////////////////

      // Do the natural
      if (p.first) {
         if (DEBUG_PRINT) std::cout << "succuss first time" << std::endl;
         max_iter = state.iters_used();
         return fn_return(state, p);
      } 
      if (DEBUG_PRINT) std::cout << "FAILURE first time" << std::endl << std::endl;


      cache_x_orig.is_l() ? fn_form_bracket_h() : fn_form_bracket_l();
      state.set_state_2();
      if (state.is_B()) {
         if (DEBUG_PRINT) std::cout << "BRACKETED PRIORITY 1" << std::endl;
         p = solve(state);
         return fn_return(state, p);
      } else {
         if (DEBUG_PRINT) std::cout << "COULD NOT BRACKET PRIORITY 1" << std::endl << std::endl;
      }

      cache_x_orig.is_l() ? fn_form_bracket_l() : fn_form_bracket_h();
      state.set_state_2();
      if (state.is_B()) {
         if (DEBUG_PRINT) std::cout << "BRACKETED PRIORITY 2" << std::endl;
         p = solve(state);
         return fn_return(state, p);
      } else {
         if (DEBUG_PRINT) std::cout << "COULD NOT BRACKET PRIORITY 2" << std::endl << std::endl;
      }

      static const char* function = "boost::math::tools::detail::solve<%3%>";
      const T x_last = state.x_next_value();
      return policies::raise_evaluation_error(function, "Unable to bracket root in boost::math::tools::detail::solve, last closest guess was %1%", x_last, boost::math::policies::policy<>());
   };

//                  from compile_test/gauss_kronrod_concept_test.cpp:11:
// ../../../libs/math/test/test_ibeta_inv.hpp(64): last checkpoint
// ../../../bin.v2/libs/math/test/test_bessel_airy_zeros.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/test_bessel_airy_zeros.o...
//  test_inverse_gaussian
//         cat "../../../bin.v2/libs/math/test/test_ibeta_inv_double.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/test_ibeta_inv_double.output"

// test_hypergeometric_laplace_transform



// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/test_vector_barycentric_rational.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/test_vector_barycentric_rational.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/cstdfloat_concept_check_1.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/compile_test/cstdfloat_concept_check_1.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/std_real_concept_check.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/compile_test/std_real_concept_check.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/std_real_concept_check_80.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/std_real_concept_check_80.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/std_real_concept_check_32.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/std_real_concept_check_32.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/std_real_concept_check_128.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/std_real_concept_check_128.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/cstdfloat_concept_check_2.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/compile_test/cstdfloat_concept_check_2.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/cstdfloat_concept_check_3.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/compile_test/cstdfloat_concept_check_3.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/cstdfloat_concept_check_4.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/compile_test/cstdfloat_concept_check_4.o...
// ...failed gcc.compile.c++ ../../../bin.v2/libs/math/test/std_real_concept_check_64.test/gcc-9/debug/link-static/threading-multi/visibility-hidden/std_real_concept_check_64.o...



// ../../../b2 --report_level=detailed --show_progress=false | grep -E "test case|failed"




// ### ??? 
//
// ### CONFIRMED TO BE REAL TESTS ###
// test_legendre
// test_ibeta_inv_double
// test_laguerre
// test_bessel_airy_zeros
// test_hermite
//
// testing.capture-output ../../../bin.v2/libs/math/example/legendre_stieltjes_example.test/gcc-9/debug/threading-multi/visibility-hidden/legendre_stieltjes_example.run






// Start here
#if 0
////////////////////////////////////////////////////////////////////////////////////////

   // Forward declare ContextF
   template <class F, class T, class Step>
   class ContextF;


   //
   //
   //
   //
   //
   template <class F, class T, class Step>
   class BoundsState {
   public:
      virtual ~BoundsState() = default;

      void set_b(BoundState<F, T, Step>& b) { b_ = &b; }

      void do_step_eval_then_x_next(ContextF<F, T, Step>& helper_f) {
         // Attempt Newton step
         const auto cx = cache_x();
         const T dx_newton = cx.calc_dx();
         const T x_next = cx.x() + dx_newton;

         // Need to do bisection
         if (!is_inbounds(x_next) || is_converge_too_slow(dx_newton)) {
            calc_next_bisection();
            return;
         }

         // Newton's step
         update_cache(helper_f(x_next));
         return;
      }

      bool is_converge_too_slow(const T& dx) const {
         return fabs(b_->dx_pp()) < fabs(dx * 4);
      }

      virtual const CacheOrder1<T, Step>& cache_x() const = 0;
      virtual NEXT<T, Step> calc_next_bisection() = 0;
      virtual void update_cache(const CacheOrder1<T, Step>& cache_x) = 0;

      bool is_inbounds(T x) const { return (l_smart() <= x) && (x <= h_smart()); }
      virtual T l_smart() const = 0;
      virtual T h_smart() const = 0;


   protected: 
      ContextF<F, T, Step>* b_;
   };



   //
   template <class F, class T, class Step>
   class ContextF {
   public:
      ContextF(F f, int digits, std::uintmax_t max_iter)
         : f_(f)
         , factor_(static_cast<T>(ldexp(1.0, 1 - digits)))
         , count_(max_iter)
         , max_iter_(max_iter) {}

      ~ContextF() {
         delete bounds_;
      }

      CacheOrder1<T, Step> operator()(T x) {
         if (0 == count_) {
            static const char* function = "boost::math::tools::detail::ContextF<%1%>";
            policies::raise_evaluation_error(function, "Unable to bracket root in boost::math::tools::detail::ContextF, last closest guess was %1%", x, boost::math::policies::policy<>());
         }
         --count_;
         return CacheOrder1<T, Step>(f_, x);
      }

      void transition_to_state(BoundsState<F, T, Step>* bounds) {
         if (this->bounds_ != nullptr)
            delete this->bounds_;
         this->bounds_ = bounds;
         this->bounds_->set_b(*this);
      }

      SolveOutput<T> solve(BoundsState<F, T, Step>* bounds) {
         transition_to_state(bounds);
         return solve();
      }

      T dx_pp() const { return dx_history_.dx_pp(); }
      int iters_used() const { return max_iter_ - count_; }

   private:
      SolveOutput<T> solve() {
         BOOST_MATH_STD_USING
         std::cout << std::scientific << std::setprecision(6) << std::endl;

         while (count_) {
            bounds_->do_step_eval_then_x_next(*this);
         }
         assert(false);
         return std::make_pair(OutputResult::failure, 0);
      }

#if 0
      std::pair<bool, T> solve(Bounds<T, Step>& bounds, CacheOrder1<T, Step> cache_x) {
         BOOST_MATH_STD_USING

         std::cout << std::scientific << std::setprecision(6) << std::endl;

         T dx_true;
         bool is_prev_newton = false;

         T x_next;
         while (count_) {
            // Calculate Newton step
            const T dx_newton = cache_x.calc_dx();
            x_next = cache_x.x() + dx_newton;  // Checked against bracket bounds below


            const bool c1 = !bounds.is_inbounds(x_next);
            if (c1 || fabs(dx_pp()) < fabs(dx_newton * 4)) {  // Bisection halves step size each iteration
               const auto p_bisect = bounds.handle_bisection(*this, cache_x);
               if (bounds.is_case_2I()) { std::cout << "is_case_2I" << std::endl;  return std::make_pair(false, 0); }

               const auto cache_x_next = p_bisect.second;
               std::cout << "STEP: BISECTION: " << " -- wanted: " << x_next << " -- got: " << cache_x_next.x() << " -- bounds: " << c1 << " -- converge: " << (fabs(dx_pp()) < fabs(dx_newton * 4)) << std::endl;
               x_next = cache_x_next.x();
               dx_history_.reset();  // Invalidate step size history

               dx_true = x_next - cache_x.x();

               cache_x = cache_x_next;
               is_prev_newton = false;
            } else {
               std::cout << "STEP: NEWTON" << std::endl;
               dx_true = x_next - cache_x.x();
               cache_x = (*this)(x_next);  // Update cache_x_ with x_next
               is_prev_newton = true;
            }

            std::cout << "AFTER STEP-- x_next: " << x_next << ", dx_true: " << dx_true << std::endl;
            bounds.print();

            assert(bounds.is_inbounds(x_next));

            // std::cout << "    UPDATED cache_x: " << std::endl;
            // std::cout << "    ";
            // cache_x.print();
            // std::cout << "    count: " << count_ << std::endl;

            // Determine if successful       
            const bool is_dx_small = fabs(cache_x.calc_dx()) < fabs(cache_x.x() * factor_);
            if (is_dx_small) {
               if (bounds.is_case_2B() || bounds.is_case_1()) {
                  std::cout << "is_case_2B && is_dx_small" << std::endl;
                  std::cout << "dx_true: " << fabs(dx_true) << ", cache_x.x(): " << fabs(cache_x.x()) << ", factor_: " << factor_ << std::endl;
                  break;  // Tolerance met
               } else {
                  bounds.print_case();
                  std::cout << "dx_true: " << fabs(dx_true) << ", cache_x.x(): " << fabs(cache_x.x()) << ", factor_: " << factor_ << std::endl;
                  std::cout << "FAILURE dx_small" << std::endl;
                  return std::make_pair(false, 0);  // Did not bracket root
               }
            }

            if (cache_x.f0() == 0) {
               std::cout << "SUCCESS f0 = 0" << std::endl;
               break;  // Function is zero
            }

            const bool is_numerical = (x_next == bounds.l_smart() || x_next == bounds.h_smart() || dx_true == 0);
            if (is_numerical) {
               if (bounds.is_case_2B()) {
                  std::cout << "is_case_2B && is_numerical" << std::endl;
                  break;  // Tolerance met
               } else {
                  std::cout << "x_next == bounds.l_smart(): " << (x_next == bounds.l_smart()) << std::endl;
                  std::cout << "x_next == bounds.h_smart(): " << (x_next == bounds.h_smart()) << std::endl;
                  std::cout << "dx_true == 0: " << (dx_true == 0) << std::endl;
                  std::cout << "FAILURE numerical" << std::endl;
                  return std::make_pair(false, 0);  // Did not bracket root
               }
            }


            if (bounds.is_case_2I()) {
               std::cout << "FAILURE is_case_2I" << std::endl;
               return std::make_pair(false, 0);
            }

            bounds.update_state(cache_x);  // Update bracket state

            if (bounds.is_case_2I()) {
               std::cout << "FAILURE is_case_2I" << std::endl;
               return std::make_pair(false, 0);
            }

            std::cout << std::endl;
            dx_history_.push(dx_newton);  // Store step size
         }

         return std::make_pair(true, cache_x.x());
      }
#endif       

   private:
      F f_;
      T factor_;
      std::uintmax_t count_;
      std::uintmax_t max_iter_;
      StepSizeHistory<T> dx_history_;
      BoundsState<F, T, Step>* bounds_;
   };



   template <class F, class T, class Step>
   class Bounds2;  // : public BoundsState<F, T, Step>;

   template <class F, class T, class Step>
   class Bounds2B;  // : public BoundsState<F, T, Step>;

   template <class F, class T, class Step>
   class Bounds2U;  // : public BoundsState<F, T, Step>;

   template <class F, class T, class Step>
   class Bounds2I;  // : public BoundsState<F, T, Step>;

   //
   //
   //
   //
   //
   template <class F, class T, class Step>
   class Bounds2 : public BoundsState<F, T, Step> {
   public:
      static Bounds2* Identify(CacheOrder1<T, Step> cache_a, CacheOrder1<T, Step> cache_b) {
         const auto cache_sorted = cache_a.sort(cache_b);
         const auto cache_l = cache_sorted.first;
         const auto cache_h = cache_sorted.second;

         if (cache_l.is_B(cache_h)) {
            return new Bounds2B<F, T, Step>(cache_l, cache_h);
         } else if (cache_l.is_U(cache_h)) {
            return new Bounds2U<F, T, Step>(cache_l, cache_h);
         } else {
            return new Bounds2I<F, T, Step>(cache_l, cache_h);
         }
      }

      NEXT<T, Step> calc_next_bisection() override {
         const auto midpoint = (cache_l_.x() + cache_h_.x()) / 2;
         const auto cache_midpoint = helper_f(midpoint);
         update_cache(cache_midpoint);
      }

      void update_cache_sign_f0(const CacheOrder1<T, Step>& cache_x) {
         if (cache_x.is_same_sign(cache_l_)) {
            cache_l_ = cache_x;
         } else {
            cache_h_ = cache_x;
         }
      }

      void update_cache_sign_dx(const CacheOrder1<T, Step>& cache_x) {
         if (cache_x.is_same_slope(cache_l_)) {  // TODO: make better
            cache_l_ = cache_x;
         } else {
            cache_h_ = cache_x;
         }
      }

      T l_smart() const override { return cache_l_.x(); }
      T h_smart() const override { return cache_h_.x(); }

      const detail::CacheOrder1<T, Step>& cache_x() const override {
         return (is_last_eval_l_) ? cache_l_ : cache_h_;
      }

      bool is_B() const { return cache_l_.is_B(cache_h_); }
      bool is_U() const { return cache_l_.is_U(cache_h_); }
      bool is_I() const { return cache_l_.is_I(cache_h_); }

      const CacheOrder1<T, Step>& cache_l() { return cache_l_; }
      const CacheOrder1<T, Step>& cache_h() { return cache_h_; }

   private:
      bool is_last_eval_l_;
      CacheOrder1<T, Step> cache_l_;
      CacheOrder1<T, Step> cache_h_;
   };

   //
   //
   //
   //
   //
   template <class F, class T, class Step>
   class Bounds2B : public Bounds2<F, T, Step> {
   public:
      Bounds2B(const CacheOrder1<T, Step>& cache_l, const CacheOrder1<T, Step>& cache_h) 
         : Bounds2<F, T, Step>(cache_l, cache_h)
      {
         assert(this->is_B());
      }

      void update_cache(const CacheOrder1<T, Step>& cache_x) override {
         update_cache_sign_f0(cache_x);
      }
   };

   //
   //
   //
   //
   //
   template <class F, class T, class Step>
   class Bounds2U : public Bounds2<F, T, Step> {
   public:
      Bounds2U(const CacheOrder1<T, Step>& cache_l, const CacheOrder1<T, Step>& cache_h) 
         : Bounds2<F, T, Step>(cache_l, cache_h)
      {
         assert(this->is_U());
      }

      void update_cache(const CacheOrder1<T, Step>& cache_x) override {
         update_cache_sign_dx(cache_x);
         if (this->is_B()) {
            this->b_->transition_to_state(
               new Bounds2B<F, T, Step>(this->cache_l_, this->cache_h_)
            );
         }
      }
   };


   //
   //
   //
   //
   //
   template <class F, class T, class Step>
   class Bounds2I : public Bounds2<F, T, Step> {
   public:
      Bounds2I(const CacheOrder1<T, Step>& cache_l, const CacheOrder1<T, Step>& cache_h) 
         : Bounds2<F, T, Step>(cache_l, cache_h)
      {
         assert(this->is_U());
      }

      void update_cache(const CacheOrder1<T, Step>& cache_x) override {
         assert(false);
      }
   };


   //
   //
   //
   //
   //
   template <class F, class T, class Step>
   class Bounds1 : public BoundsState<F, T, Step> {
   public:
      Bounds1(T l, const detail::CacheOrder1<T, Step>& cache_x, T h) 
         : cache_x_(cache_x)
      {
         limit_ = cache_x_.is_l() ? h : l;
      }

      NEXT<T, Step> calc_next_bisection() override {
         const auto cache_limit = (*this->b_)(limit_);
         auto* bounds_new = Bounds2<F, T, Step>::Identify(cache_limit, cache_x_);
         this->b_->transition_to_state(bounds_new);
      }

      T l_smart() const override { return cache_x_.is_l() ? cache_x_.x() : limit_; }
      T h_smart() const override { return cache_x_.is_h() ? cache_x_.x() : limit_; }

      void update_cache(const CacheOrder1<T, Step>& cache_x) override { cache_x_ = cache_x; }

      const CacheOrder1<T, Step>& cache_x() const override { return cache_x_; }

   private:
      CacheOrder1<T, Step> cache_x_;
      T limit_;
   };



///////////////////////////////////////////////////////////////////////


   template <class F, class T, class Step>
   T solve(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>()))) {
      // Create ContextF
      ContextF<F, T, Step> b(f, digits, max_iter);

      // Create cache_x
      const auto cache_x = b(x);
      if (0 == cache_x.f0()) return x;  // Check if root is already found
      
      // Create Bounds
      auto bounds = Bounds1<F, T, Step>(xl, cache_x, xh);

      // Solve 
      const auto p = b.solve(&bounds);

      return p.second;

#if 0
      if (p.first) {
         std::cout << "succuss first time" << std::endl;
         max_iter = b.iters_used();
         return p.second;
      } else {
         std::cout << "FAILURE first time" << std::endl << std::endl;

         bounds = Bounds<T, Step>(xl, cache_x, xh);
         cache_x.is_l() ? bounds.eval_h_orig(b) : bounds.eval_l_orig(b);
         if (bounds.is_case_2B()) { 
            std::cout << "BRACKETED PRIORITY 1" << std::endl;
            return b.solve(bounds, cache_x).second;
         } else {
            std::cout << "COULD NOT BRACKET PRIORITY 1" << std::endl << std::endl;
         }

         bounds = Bounds<T, Step>(xl, cache_x, xh);
         cache_x.is_l() ? bounds.eval_l_orig(b) : bounds.eval_h_orig(b);
         if (bounds.is_case_2B()) { 
            std::cout << "BRACKETED PRIORITY 2" << std::endl;
            return b.solve(bounds, cache_x).second;
         } else {
            std::cout << "l:    " << bounds.l_smart()      << ", h:    " << bounds.h_smart() << std::endl;
            std::cout << "f0_l: " << bounds.cache_l().f0() << ", f0_h: " << bounds.cache_h().f0() << std::endl;
            std::cout << "COULD NOT BRACKET PRIORITY 2" << std::endl << std::endl;
         }
      }
#endif

      static const char* function = "boost::math::tools::detail::ContextF<%1%>";
      policies::raise_evaluation_error(function, "Unable to bracket root in boost::math::tools::detail::ContextF, last closest guess was %1%", x, boost::math::policies::policy<>());
      return 0;
   };






#if 0


   //
   //
   //
   //
   //
   template <class T, class Step>
   class Bounds2 {
   public:
      Bounds2(T l, const detail::CacheOrder1<T, Step>& cache_x, T h) 
         : var_l_(LimitOrig<T>(l)), var_h_(LimitOrig<T>(h))
      {
         if (cache_x.is_l()) { 
            is_last_eval_l_ = true;
            var_l_ = cache_x;
         } else {
            is_last_eval_l_ = false;
            var_h_ = cache_x;
         }
      }

      template <class H>
      void do_step_eval_then_x_next(H& helper_f) {
         const auto cx = cache_x();
         const T x_next = cx.x() + cx.calc_dx();

         if (!is_inbounds(x_next) || is_converge_too_slow()) {
            return calc_next_bisection(helper_f);
         }

         // Newton's step
         update_cache(helper_f(x_next));
      }

      template <class H>
      NEXT<T, Step> calc_next_bisection(H& helper_f) {
         if (is_case_1()) {
            materialize(helper_f);
            return;
         }
         const T x_mid = (x_l() + x_h()) / 2;
         update_cache_2(helper_f(x_mid));
         return;
      }

      void update_cache(const detail::CacheOrder1<T, Step>& cache_x) {
         if (is_case_1()) {
            update_cache_1(cache_x);
         } else {
            update_cache_2(cache_x);
         }
      }

      void update_cache_1(const detail::CacheOrder1<T, Step>& cache_x) {
         update_cache_using_dx(cache_x);
         if (size() == 2) {
            determine_state_2();
         }
      }

      void update_cache_2(const detail::CacheOrder1<T, Step>& cache_x) {
         if (is_case_2B()) {
            update_cache_using_sign(cache_x);
         } else {
            assert(is_case_2U());
            update_cache_2U(cache_x);
         }
      }

      void update_cache_2U(const detail::CacheOrder1<T, Step>& cache_x) {
         const T dx_l_orig = cache_l().calc_dx();
         const T dx_h_orig = cache_h().calc_dx();

         update_cache_using_dx(cache_x);

         if (cache_l().is_B(cache_h())) {  // Now bracket root!
            state_ = BracketState::case_2B;
         } else {
            const T dx_l = cache_l().calc_dx();
            const T dx_h = cache_h().calc_dx();

            if (dx_l_orig < dx_l || dx_h_orig < dx_h) {
               state_ = BracketState::case_2I;
            }
         }
      }

      template <class H>
      void materialize(H& helper_f) {
         // If var_l_ holds a LimitOrig, then we need to evaluate it
         if (var_l_.is_limit_orig()) {
            set_l(helper_f(var_l_.x()));
         }
         // If var_h_ holds a LimitOrig, then we need to evaluate it
         if (var_h_.is_limit_orig()) {
            set_h(helper_f(var_h_.x()));
         }
         determine_state_2();
      }

      void determine_state_2() {
         assert(size() == 2);
         if (cache_l().is_B(cache_h())) {
            state_ = BracketState::case_2B;
         } else if (cache_l().is_U(cache_h())) {
            state_ = BracketState::case_2U;
         } else {
            state_ = BracketState::case_2I;
         }
      }

      void update_state_using_dx(const CacheOrder1<T, Step>& cache_x) {
         if (cache_x.is_l()) {
            var_l_ = cache_x;
         } else {
            var_h_ = cache_x;
         }
      }

      void update_state_using_sign(const CacheOrder1<T, Step>& cache_x) {
         if (has_l()) {
            if (cache_x.is_same_sign(cache_l())) {
               var_l_ = cache_x;
            } else {
               var_h_ = cache_x;
            }
         } else {
            if (cache_x.is_same_sign(cache_h())) {
               var_h_ = cache_x;
            } else {
               var_l_ = cache_x;
            }
         }
      }

////////////////////////////////////////////////////////////////////////////////

      detail::CacheOrder1<T, Step> cache_x() const { return is_last_eval_l_ ? val_l_ : var_h_; }

      int size() const { return var_l_.is_cache() + var_h_.is_cache(); }

      bool is_inbounds(T x) const { return (l() <= x) && (x <= h()); }

      void set_l(const detail::CacheOrder1<T, Step>& cache_x) {
         is_last_eval_l_ = true;
         var_l_ = cache_x;
      }
      void set_h(const detail::CacheOrder1<T, Step>& cache_x) {
         is_last_eval_l_ = false;
         var_h_ = cache_x;
      }

      T l() const { return var_l_.x(); }
      T h() const { return var_h_.x(); }

   private:
      BracketState state_;
      bool is_last_eval_l_;
      booth::variant<LimitOrig<T>, detail::CacheOrder1<T, Step>> var_l_;
      booth::variant<LimitOrig<T>, detail::CacheOrder1<T, Step>> var_h_;
      StepSizeHistory<T> dx_history_;
   };


















   //
   template <class T, class Step>
   class Bounds {
   public:
      Bounds(T l, detail::CacheOrder1<T, Step> cache_x, T h) 
         : state_(BracketState::case_1)
         , l_orig_(l)
         , h_orig_(h)
      {
         update_state_using_dx(cache_x);
      }

      void update_state_using_dx(const CacheOrder1<T, Step>& cache_x) {
         if (cache_x.is_l()) {
            opt_cache_l_ = cache_x;
         } else {
            opt_cache_h_ = cache_x;
         }
      }

      void update_state_using_sign(const CacheOrder1<T, Step>& cache_x) {
         if (has_l()) {
            if (cache_x.is_same_sign(cache_l())) {
               opt_cache_l_ = cache_x;
            } else {
               opt_cache_h_ = cache_x;
            }
         } else {
            if (cache_x.is_same_sign(cache_h())) {
               opt_cache_h_ = cache_x;
            } else {
               opt_cache_l_ = cache_x;
            }
         }
      }

      void determine_state_exhaustive() {
         if (cache_l().is_B(cache_h())) {
            state_ = BracketState::case_2B;
         } else if (cache_l().is_U(cache_h())) {
            state_ = BracketState::case_2U;
         } else {
            state_ = BracketState::case_2I;
         }
      }

      bool update_state(const CacheOrder1<T, Step>& cache_x) {
         if (is_case_1()) {
            const bool c1 = cache_x.is_l() && has_h();
            const bool c2 = cache_x.is_h() && has_l();
            update_state_using_dx(cache_x);
            if (c1 || c2) {
               determine_state_exhaustive();
            }
         } else if (is_case_2B()) {
            update_state_using_sign(cache_x);
            if (!cache_l().is_B(cache_h())) {
               std::cout << std::endl;
               std::cout << "SOMETHING IS WRONG WRONG WRONG" << std::endl;
               print();
               cache_x.print();
               assert(false);
            }
         } else {
            assert(is_case_2U());
            update_state_2U(cache_x);
         } 
         return true;
      }

      void update_state_2U(const CacheOrder1<T, Step>& cache_x) {
         const T dx_l_orig = cache_l().calc_dx();
         const T dx_h_orig = cache_h().calc_dx();

         update_state_using_dx(cache_x);

         if (cache_l().is_B(cache_h())) {  // Now bracket root!
            state_ = BracketState::case_2B;
         } else {
            const T dx_l = cache_l().calc_dx();
            const T dx_h = cache_h().calc_dx();

            if (dx_l_orig < dx_l || dx_h_orig < dx_h) {
               state_ = BracketState::case_2I;
            }
         }
      }

      static T midpoint_smart(T l, T h) {
         BOOST_MATH_STD_USING

         const T t_min = std::numeric_limits<T>::min();
         const T t_max = std::numeric_limits<T>::max();
         const T t_lowest = std::numeric_limits<T>::lowest();
         const T t_inf = std::numeric_limits<T>::infinity();

         const T l_abs = fabs(l);
         const T h_abs = fabs(h);
         if (std::max(l_abs, h_abs) < 100 || // Numbers are small
            (h - l) < std::min(l_abs, h_abs)) { // Numbers are close
            return (l + h) / 2;
         }

         // Want both numbers to have the same sign
         if (sign(l) * sign(h) == -1) { return 0; }

         // If l is -Inf make it the lowest float
         if (l == -t_inf) {
            l = t_lowest;
         } else if (h == t_inf) {
            h = t_max;
         }

         // If l is 0 make it the smallest float
         if (l == 0) {
            l = t_min;
         } else if (h == 0) {
            h = -t_min, l;
         }

         return copysign(sqrt(l * h), l);
      }

      template <class H>
       std::pair<bool, CacheOrder1<T, Step>> handle_bisection(H& helper_f, const CacheOrder1<T, Step>& cache_x) {
         if (is_case_1()) {
            const T x_next = cache_x.is_l() ? h_orig_ : l_orig_;
            
            CacheOrder1<T, Step> cache_x_next;

            if (cache_x.is_l()) {
               cache_x_next = helper_f(x_next);
               opt_cache_h_ = cache_x_next;
            } else {
               cache_x_next = helper_f(x_next);
               opt_cache_l_ = cache_x_next;
            }

            determine_state_exhaustive();
            print_case();
         }

         // const T x_mid = (l_smart() + h_smart()) / 2;
         const T x_mid = midpoint_smart(l_smart(), h_smart());

         std::cout << "x_mid: " << x_mid << std::endl;

         const auto cache_x_next = helper_f(x_mid);

         std::cout << "post eval" << std::endl;

         return std::make_pair(true, cache_x_next);
      }

      template <class H>
      void eval_l_orig(H& helper_f) {
         const auto cache_x = helper_f(l_orig_);
         push_front(cache_x);
         determine_state_exhaustive();
      }
      template <class H>
      void eval_h_orig(H& helper_f) {
         const auto cache_x = helper_f(h_orig_);
         push_back(cache_x);
         determine_state_exhaustive();
      }

      void push_front(const CacheOrder1<T, Step>& cache_x) {
         assert(size() == 1);
         if (has_l()) { opt_cache_h_ = opt_cache_l_; }
         opt_cache_l_ = cache_x;
      }
      void push_back(const CacheOrder1<T, Step>& cache_x) {
         assert(size() == 1);
         if (has_h()) { opt_cache_l_ = opt_cache_h_; }
         opt_cache_h_ = cache_x;
      }

      int size() const { return has_l() + has_h(); }
      bool has_l() const { return opt_cache_l_.has_data(); }
      bool has_h() const { return opt_cache_h_.has_data(); }

      const CacheOrder1<T, Step>& cache_l() const { return opt_cache_l_.c(); }
      const CacheOrder1<T, Step>& cache_h() const { return opt_cache_h_.c(); }

      T l_smart() const { return has_l() ? cache_l().x() : l_orig_; }
      T h_smart() const { return has_h() ? cache_h().x() : h_orig_; }

      bool is_inbounds(T x) const { return (l_smart() <= x) && (x <= h_smart()); }

      bool is_case_1() const { return state_ == BracketState::case_1; }
      bool is_case_2U() const { return state_ == BracketState::case_2U; }
      bool is_case_2B() const { return state_ == BracketState::case_2B; }
      bool is_case_2I() const { return state_ == BracketState::case_2I; }

      std::string get_case_string() const {
         if (is_case_1()) { return "case_1"; }
         if (is_case_2U()) { return "case_2U"; }
         if (is_case_2B()) { return "case_2B"; }
         assert(is_case_2I());
         return "case_2I";
      }

      void print_case() const { std::cout << get_case_string() << std::endl; }

      void print() const {
         std::cout << "BOUNDS: " << get_case_string() << std::endl;
         std::cout << "  l_orig: " << l_orig_ << std::endl;
         if (has_l()) { cache_l().print(); }
         if (has_h()) { cache_h().print(); }
         std::cout << "  h_orig: " << h_orig_ << std::endl;
      }

   private:
      BracketState state_;
      T l_orig_;
      detail::Optional<detail::CacheOrder1<T, Step>> opt_cache_l_;
      detail::Optional<detail::CacheOrder1<T, Step>> opt_cache_h_;
      T h_orig_;
   };

   template <class F, class T, class Step>
   T solve(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>()))) {
      // Create ContextF
      ContextF<F, T, Step> helper_f(f, digits, max_iter);

      // Create cache_x
      const auto cache_x = helper_f(x);
      if (0 == cache_x.f0()) return x;  // Check if root is already found
      
      // Create Bounds2
      auto bounds = Bounds1<T, Step>(xl, cache_x, xh);

      // Do Newton loop
      const auto p = helper_f.solve(bounds, cache_x);

      if (p.first) {
         std::cout << "succuss first time" << std::endl;
         max_iter = helper_f.iters_used();
         return p.second;
      } else {
         std::cout << "FAILURE first time" << std::endl << std::endl;

         bounds = Bounds2<T, Step>(xl, cache_x, xh);
         cache_x.is_l() ? bounds.eval_h_orig(helper_f) : bounds.eval_l_orig(helper_f);
         if (bounds.is_case_2B()) { 
            std::cout << "BRACKETED PRIORITY 1" << std::endl;
            return helper_f.solve(bounds, cache_x).second;
         } else {
            std::cout << "COULD NOT BRACKET PRIORITY 1" << std::endl << std::endl;
         }

         bounds = Bounds2<T, Step>(xl, cache_x, xh);
         cache_x.is_l() ? bounds.eval_l_orig(helper_f) : bounds.eval_h_orig(helper_f);
         if (bounds.is_case_2B()) { 
            std::cout << "BRACKETED PRIORITY 2" << std::endl;
            return helper_f.solve(bounds, cache_x).second;
         } else {
            std::cout << "l:    " << bounds.l_smart()      << ", h:    " << bounds.h_smart() << std::endl;
            std::cout << "f0_l: " << bounds.cache_l().f0() << ", f0_h: " << bounds.cache_h().f0() << std::endl;
            std::cout << "COULD NOT BRACKET PRIORITY 2" << std::endl << std::endl;
         }
      }

      static const char* function = "boost::math::tools::detail::ContextF<%1%>";
      policies::raise_evaluation_error(function, "Unable to bracket root in boost::math::tools::detail::ContextF, last closest guess was %1%", x, boost::math::policies::policy<>());
      return 0;
   };

#if 0
   template <class F, class T, class Step>
   T solve(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>()))) {
      // Create ContextF
      ContextF<F, T, Step> helper_f(f, digits, max_iter);

      // Create cache_x
      const auto cache_x = helper_f(x);
      if (0 == cache_x.f0()) return x;  // Check if root is already found
      
      // Create Bounds
      auto bounds = Bounds<T, Step>(xl, cache_x, xh);

      // Do Newton loop
      const auto p = helper_f.solve(bounds, cache_x);

      if (p.first) {
         std::cout << "succuss first time" << std::endl;
         max_iter = helper_f.iters_used();
         return p.second;
      } else {
         std::cout << "FAILURE first time" << std::endl << std::endl;

         bounds = Bounds<T, Step>(xl, cache_x, xh);
         cache_x.is_l() ? bounds.eval_h_orig(helper_f) : bounds.eval_l_orig(helper_f);
         if (bounds.is_case_2B()) { 
            std::cout << "BRACKETED PRIORITY 1" << std::endl;
            return helper_f.solve(bounds, cache_x).second;
         } else {
            std::cout << "COULD NOT BRACKET PRIORITY 1" << std::endl << std::endl;
         }

         bounds = Bounds<T, Step>(xl, cache_x, xh);
         cache_x.is_l() ? bounds.eval_l_orig(helper_f) : bounds.eval_h_orig(helper_f);
         if (bounds.is_case_2B()) { 
            std::cout << "BRACKETED PRIORITY 2" << std::endl;
            return helper_f.solve(bounds, cache_x).second;
         } else {
            std::cout << "l:    " << bounds.l_smart()      << ", h:    " << bounds.h_smart() << std::endl;
            std::cout << "f0_l: " << bounds.cache_l().f0() << ", f0_h: " << bounds.cache_h().f0() << std::endl;
            std::cout << "COULD NOT BRACKET PRIORITY 2" << std::endl << std::endl;
         }
      }

      static const char* function = "boost::math::tools::detail::ContextF<%1%>";
      policies::raise_evaluation_error(function, "Unable to bracket root in boost::math::tools::detail::ContextF, last closest guess was %1%", x, boost::math::policies::policy<>());
      return 0;
   };
#endif
#endif
#endif

}  // namespace detail




template <class F, class T>
T newton_raphson_iterate(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   return detail::solve<F, T, detail::StepNewton<T>>(f, x, xl, xh, digits, max_iter);
}

template <class F, class T>
inline T newton_raphson_iterate(F f, T guess, T min, T max, int digits) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::uintmax_t m = (std::numeric_limits<std::uintmax_t>::max)();
   return newton_raphson_iterate(f, guess, min, max, digits, m);
}



#if 1

template <class F, class T>
T newton_raphson_iterate_classic(F f, T guess, T min, T max, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   BOOST_MATH_STD_USING

   static const char* function = "boost::math::tools::newton_raphson_iterate_classic<%1%>";
   if (min > max)
   {
      return policies::raise_evaluation_error(function, "Range arguments in wrong order in boost::math::tools::newton_raphson_iterate_classic(first arg=%1%)", min, boost::math::policies::policy<>());
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
   std::cout << "newton_raphson_iterate_classic, guess = " << guess << ", min = " << min << ", max = " << max
      << ", digits = " << digits << ", max_iter = " << max_iter << "\n";
#endif

   do {
      last_f0 = f0;
      delta2 = delta1;
      delta1 = delta;
      std::cout << std::endl;
      std::cout << "EVAL IS: " << count << ", at: " << result << std::endl;
      detail::unpack_tuple(f(result), f0, f1);
      --count;

      std::cout << "result: " << result << ", f0: " << f0 << ", f1: " << f1;

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
         std::cout << ", delta: " << delta <<  std::endl;
      }
#ifdef BOOST_MATH_INSTRUMENT
      std::cout << "Newton iteration " << max_iter - count << ", delta = " << delta << ", residual = " << f0 << "\n";
#endif
      if (fabs(delta * 2) > fabs(delta2))
      {
         // Last two steps haven't converged.
         T shift = (delta > 0) ? (result - min) / 2 : (result - max) / 2;
         if ((result != 0) && (fabs(shift) > fabs(result)))
         {
            delta = sign(delta) * fabs(result) * 1.1f; // Protect against huge jumps!
            //delta = sign(delta) * result; // Protect against huge jumps! Failed for negative result. https://github.com/boostorg/math/issues/216
         }
         else
            delta = shift;
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
         return policies::raise_evaluation_error(function, "There appears to be no root to be found in boost::math::tools::newton_raphson_iterate_classic, perhaps we have a local minima near current best guess of %1%", guess, boost::math::policies::policy<>());
      }
   }while(count && (fabs(result * factor) < fabs(delta)));

   std::cout << "delta:               " << delta << std::endl;
   std::cout << "fabs(result * factor)" << fabs(result * factor) << std::endl;
   std::cout << "fabs(delta)          " << fabs(delta) << std::endl;

   max_iter -= count;

#ifdef BOOST_MATH_INSTRUMENT
   std::cout << "Newton Raphson required " << max_iter << " iterations\n";
#endif

   return result;
}

template <class F, class T>
inline T newton_raphson_iterate_classic(F f, T guess, T min, T max, int digits) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::uintmax_t m = (std::numeric_limits<std::uintmax_t>::max)();
   return newton_raphson_iterate_classic(f, guess, min, max, digits, m);
}


#endif













#else 



/////////////////////////////////


   // Enforces bracketing when root finding with Newton's method
   template <class F, class T, class Step>
   class BracketHelper {
      // Stores the last two Newton step sizes
      class StepSizeHistory {
      public:
         StepSizeHistory() { reset(); }

         void reset() {
            dx_p_ = (std::numeric_limits<T>::max)();
            dx_pp_ = (std::numeric_limits<T>::max)();
         }

         void push(T input) {
            dx_pp_ = dx_p_;
            dx_p_ = input;
         }

         T dx_p() const { return dx_p_; }
         T dx_pp() const { return dx_pp_; }

      private:
         T dx_p_;
         T dx_pp_;
      };

      // 
      

      

      static Bounds create_bracket_candidate(ContextF<F, T>& helper_f, const detail::CacheOrder1<T>& cache_x, T bound) {
         if (is_can_eval_bound(bound)) {
            return Bounds(cache_x, helper_f(bound));
         } else {
            return rollout_for_bracket_candidate(helper_f, cache_x, bound);
         }
      }

      static

   public:
      static T solve(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>()))) {
         // // Create ContextF
         // ContextF<F, T> helper_f(f, digits, max_iter);

         // // Create cache_x
         // const auto cache_x = helper_f(x);

         // // Check if root is already found
         // if (0 == cache_x.f0()) return x;

         // // Get preferred bound
         // const T dx = Step().calc_dx(cache_x);
         // const T bound_pref = (0 <= dx) ? xl : xh;

         // // Check if bound is valid
         // const auto bc_pref = create_bracket_candidate(helper_f, cache_x, bound_pref);




         const bool is_can_eval_bound = is_can_eval_bound(bound);
         if (is_can_eval_bound) {
            const auto cache_b = detail::CacheOrder1<T>(f, bound);
            const auto bc = Bounds(cache_x, cache_b);
            if (bc.may_contain_root()) {
               return solve_bc(f, bc, digits, max_iter, cache_x);
            } else {
               return std::make_pair(false, 0);
            }
         }


         if (0 <= dx) {  // Want to go lower
            const auto p_l = solve_l(f, x, xl, xh, digits, max_iter, cache_x);
         } else {

         }
      }

      static std::pain<bool,T> solve_bound_prefered(F f, BranchCandidate bc, int digits, std::uintmax_t& max_iter) {
         if (bc.can_bracket()) {
            return solve_bracket(f, bc, digits, max_iter, cache_x);
         } else if (bc.can_brent()) {
            return solve_brent(f, bc, digits, max_iter, cache_x);
         } else {
            return std::make_pair(false, 0);
         }
      }

      static std::pair<bool,T> solve_l(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter, const detail::CacheOrder1<T>& cache_x) {
         const bool is_can_eval_bound(xl);
         if (is_can_eval_bound) {
            const auto bc = BranchCandidate(cache_x, detail::CacheOrder1<T>(f, xl));
            if (bc.may_contain_root()) {
               return solve_bc(f, x, xl, xh, digits, max_iter, cache_x);
            } else {
               return std::make_pair(false, 0);
            }
         }

         if (!is_can_eval_bound) {
            // Find lower bound
         }


      }




         // std::cout << "================= solve ====================" << std::endl;

         // const auto cache_debug = detail::CacheOrder1<T>(f, x);
         // std::cout << "x: " << cache_debug.x() << ", f0: " << cache_debug.f0() << ", f1: " << cache_debug.f1() << std::endl;
         // // x: 0.56763, f0: -0.00218081, f1: 0.0721978

         // const auto cache_debug_l = detail::CacheOrder1<T>(f, xl);
         // std::cout << "xl: " << cache_debug_l.x() << ", f0: " << cache_debug_l.f0() << ", f1: " << cache_debug_l.f1() << std::endl;
         // // xl: 0, f0: -0.999001, f1: 0

         // const auto cache_debug_y = detail::CacheOrder1<T>(f, 1.0);
         // std::cout << "y: " << cache_debug_y.x() << ", f0: " << cache_debug_y.f0() << ", f1: " << cache_debug_y.f1() << std::endl;
         // // y: 1, f0: 0.000998974, f1: 1.64892e-07

         // const auto cache_debug_z = detail::CacheOrder1<T>(f, 2.0);
         // std::cout << "z: " << cache_debug_z.x() << ", f0: " << cache_debug_z.f0() << ", f1: " << cache_debug_z.f1() << std::endl;
         // // z: 2, f0: 0.000998974,flip f1: 2.88776e-33


         // std::cout << "xh: " << xh << std::endl;
         // // xh: 3.40282e+38


         // Get maximum value of T
         // const T x_limit = std::numeric_limits<T>::max();

         // std::cout << "x_limit: " << x_limit << std::endl;
         // const bool bb = x_limit == xh;
         // std::cout << bb << std::endl;

         // assert(false);

         const auto cache_debug_h = detail::CacheOrder1<T>(f, xh);
         std::cout << "xh: " << cache_debug_h.x() << ", f0: " << cache_debug_h.f0() << ", f1: " << cache_debug_h.f1() << std::endl;
         // Did not finish


         BracketHelper helper(f, x, xl, xh, digits, max_iter);

         T x_sol = helper.solve_impl();
         max_iter -= helper.count();  // Make max_iter the number of iterations used.
         std::cout << "max_iter: " << max_iter << std::endl;
         return x_sol;
      }

   private:
      // Same as the inputs to newton_raphson_iterate
      BracketHelper(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter)
         : f_(f)
         , cache_x_(f, x)
         , cache_xl_(f, xl)
         , cache_xh_(f, xh)
         , count_(max_iter)
         , factor_(static_cast<T>(ldexp(1.0, 1 - digits)))  // factor_ = 1.0 * 2^(1-digits)
      {
         if (xh < xl) {
            static const char* function = "boost::math::tools::detail::BracketHelper<%1%>";
            policies::raise_evaluation_error(function, "Range arguments in wrong order in boost::math::tools::detail::BracketHelper(first arg=%1%)", xl, boost::math::policies::policy<>());
         }

         std::cout << "================= factor_: " << factor_ << "====================" << std::endl;
      }


      T calc_dx_next_default() const { return Step().calc_dx(cache_x_); }
      T calc_x_next_default() const { return x() - calc_dx_next_default(); }

      bool can_bracket_l() const { return f0() * f0_l() <= 0; }
      bool can_bracket_h() const { return f0() * f0_h() <= 0; }

      T solve_bracket_l() {
         set_cache_H_to_cache_x();  // Evaluated point x becomes new high bound
         return solve_bracket();
      }
      T solve_bracket_h() {
         set_cache_L_to_cache_x();  // Evaluated point x becomes new low bound
         return solve_bracket();
      }

      bool want_to_go_l() const { return 0 <= calc_dx_next_default(); }  // Because of negative sign in update formula x_next = x() - dx

      //
      //
      //
      //
      //
      T solve_impl() {
         if (0 == f0()) return x();  // Solved at initial x


      }


#if 0
      T solve_impl() {
         if (0 == f0()) return x();  // Solved at initial x

         // Calculate direction towards root
         const T dx = Step().calc_dx(cache_x_);
         const T x_next = x() - dx;  // x_next calculated with step dx

         const bool can_bracket_l = f0() * f0_l() <= 0;  // Sign flip between x and xl
         const bool can_bracket_h = f0() * f0_h() <= 0;  // Sign flip between x and xh
         const bool want_to_go_l = 0 <= dx;  // Because of negative sign in update formula x_next = x() - dx

         if (can_bracket_l && (want_to_go_l || !can_bracket_h)) {
            set_cache_H_to_cache_x();  // Evaluated point x becomes new high bound
            return solve_bracket(x_next);
         } else if (can_bracket_h) {
            set_cache_L_to_cache_x();  // Evaluated point x becomes new low bound
            return solve_bracket(x_next);
         } else {
            static const char* function = "boost::math::tools::detail::BracketHelper<%1%>";
            return policies::raise_evaluation_error(function, "There appears to be no root to be found in boost::math::tools::detail::BracketHelper, perhaps we have a local minima near current best guess of %1%", x(), boost::math::policies::policy<>());
         }
      }
#else

////////////////////////////////////

      T solve_impl() {
         if (0 == f0()) return x();  // Solved at initial x

         std::cout << "SOLVE_IMPL" << std::endl;

         // Bracket root in direction of dx if possible
         if (want_to_go_l() && can_bracket_l()) {
            std::cout << "DIRECT --> SOLVE_BRACKET_L" << std::endl;
            return solve_bracket_l();
         } else if (can_bracket_h()) {
            std::cout << "DIRECT --> SOLVE_BRACKET_H" << std::endl;
            return solve_bracket_h();
         }

         return solve_classic_with_fallback();
      }

      T solve_classic_with_fallback() {
         std::cout << "Doing classic with fallback" << std::endl;

         // Copy caches
         const auto cache_x_copy = cache_x_;
         const auto cache_xh_copy = cache_xh_;
         const auto cache_xl_copy = cache_xl_;

         bool is_success_classic = false;
         const T root = solve_classic(&is_success_classic);

         if (is_success_classic) {
            return root;
         } else {
            // Restore caches
            cache_x_ = cache_x_copy;
            cache_xh_ = cache_xh_copy;
            cache_xl_ = cache_xl_copy;

            // Try bracketing
            if (can_bracket_l()) {
               return solve_bracket_l();
            } else if (can_bracket_h()) {
               return solve_bracket_h();
            }
         }

         static const char* function = "boost::math::tools::detail::BracketHelper<%1%>";
         return policies::raise_evaluation_error(function, "There appears to be no root to be found in boost::math::tools::detail::BracketHelper, perhaps we have a local minima near current best guess of %1%", x(), boost::math::policies::policy<>());
      }

      //
      // The classic solver is used when the root is not bracketed. If the classic solver detects
      // a sign flip it passes control to the bracket solver.
      //
      T solve_classic(bool* is_success) {
         BOOST_MATH_STD_USING
         static const char* function = "boost::math::tools::detail::BracketHelper<%1%>";

         const bool same_sign_lh = 1 == sign(f0_l()) * sign(f0_h());
         const bool same_sign_lx = 1 == sign(f0_l()) * sign(f0());

         if (!same_sign_lh || !same_sign_lx) {
            *is_success = false;
            // return policies::raise_evaluation_error(function, "There appears to be no root to be found in 
            // boost::math::tools::newton_raphson_iterate, perhaps we have a local minima near current best guess 
            // of %1%", guess, boost::math::policies::policy<>());
            return policies::raise_evaluation_error(function, "Not same sign boost::math::tools::detail::BracketHelper %1%", xl(), boost::math::policies::policy<>());
         }

         const bool is_negate = (1 == sign(f0_l())) ? false : true;

         const auto fn_unoriented = [&](T x) {
            T brent_0;
            T brent_1;
            detail::unpack_tuple(f_(x), brent_0, brent_1);
            std::cout << " BRENT eval at: " << x << ", f(x): " << brent_0 << std::endl;
            return brent_0;
         };

         const auto fn_brent = [&](T x) {
            T brent_0 = fn_unoriented(x);
            return is_negate ? -brent_0 : brent_0;
         };


         for (double i = 0.0; i <= 1.0; i += 0.0001) {
            // std::cout << "fn_brent(" << i << "): ";  //  << fn_brent(i) << std::endl;
            fn_brent(i);
         }

         std::cout << std::endl;


         const auto pair_x_fx = brent_find_minima(fn_brent, xl(), xh(), 1024);
         const auto x_min = pair_x_fx.first;
         const auto fx_min = fn_unoriented(x_min);

         std::cout << "x_min: " << x_min << std::endl;
         std::cout << "fx_min: " << fx_min << std::endl;

         if (1 == sign(f0_l()) * sign(fx_min)) {
            std::cout << "CLASSIC could not find root bracket" << std::endl;
            std::cout << "x:  " << x()  << ", f0:   " << f0()   << ", f1: " << cache_x_.f1() << std::endl;
            std::cout << "xl: " << xl() << ", f0_l: " << f0_l() << ", f1_l: " << f1_l() << std::endl;
            std::cout << "xh: " << xh() << ", f0_h: " << f0_h() << ", f1_h: " << f1_h() << std::endl;
            *is_success = false;
            return x_min;  // Return doesn't matter
         } else {
            *is_success = true;
            std::cout << "CLASSIC found root bracket" << std::endl;
            cache_x_ = detail::CacheOrder1<T>(f_, x_min);  // Update cache_x_ with x_next



            return (sign(fx_min) == 1) ? solve_bracket_h() : solve_bracket_l();
         }
      }

#endif

      T solve_bracket() {
         BOOST_MATH_STD_USING
         static const char* function = "boost::math::tools::detail::BracketHelper<%1%>";

         std::cout << "SOLVE_BRACKET" 
                   << "RealType: " << typeid(T).name() 
                   << " -- count: " << count_ << "====================" << std::endl;


         if (f0_l() == 0.0) {  // Low bound is root
            std::cout << "Low bound is root" << std::endl;
            return xl();
         } else if (f0_h() == 0.0) {  // High bound is root
            std::cout << "High bound is root" << std::endl;
            return xh();
         } else if (0 < f0_l() * f0_h()) { // Bounds do not bracket root
            std::cout << "Bounds do not bracket root" << std::endl;
            return policies::raise_evaluation_error(function, "Does not bracket boost::math::tools::detail::BracketHelper %1%", xl(), boost::math::policies::policy<>());
         }
         
         // Orient cache_l and cache_h so that f0_l is negative and f0_h is positive
         if (0 < f0_l()) std::swap(cache_xl_, cache_xh_);

         T x_next;
         do {
            count_--;

            // Calculate Newton step
            const T dx_newton = Step().calc_dx(cache_x_);
            x_next = x() - dx_newton;  // Checked against bracket bounds below

            if (!is_x_between_l_and_h(x_next) ||
               fabs(dx_pp()) < fabs(dx_newton * 4)) {  // Bisection halves step size each iteration

               // if xl and xh have different signs set x_next to 0
               if (-1 == sign(xl()) * sign(xh())) {
                  x_next = 0;
               } else {
                  const T x_S = (fabs(xl()) < fabs(xh())) ? xl() : xh();
                  const T x_B = (fabs(xl()) < fabs(xh())) ? xh() : xl();
                  const T SB_mag = fabs(xl() - xh());

                  if (x_S < SB_mag && 1 < fabs(x_S)) {  // Spans orders of magnitude
                     x_next = copysign(sqrt(fabs(x_S)) * sqrt(fabs(x_B)), x_S);
                  } else {
                     x_next = (xl() + xh()) / 2;  // Midpoint
                  }
               }

               dx_history_.reset();  // Invalidate step size history
            }

            const T dx = x_next - x();
            cache_x_ = detail::CacheOrder1<T>(f_, x_next);  // Update cache_x_ with x_next

            std::cout << "x: " << x() << ", dx: " << dx << ", dx_newton: " << dx_newton 
                      << ", f1: " << cache_x_.f1() << ", f0: " << cache_x_.f0()
                      << ", xl: " << xl() << ", xh: " << xh() 
                      << std::endl;

            if (fabs(dx) < fabs(x() * factor_)) {
               std::cout << "EXIT: tolerance met" << std::endl;
               break; }  // Tolerance met
            if (x_next == xl() || x_next == xh() || dx == 0) { 
               std::cout << "EXIT: numerical reason" << std::endl;
               break; }
               
            
            // Update bracket
            (f0() < 0.0) ? set_cache_L_to_cache_x() : set_cache_H_to_cache_x();

            dx_history_.push(dx);  // Store step size

            if (count_ == 0) {
               std::cout << "EXIT: max_iter reached" << std::endl;
               break;
            }
         } while (count_);

         return x();
      }

      bool is_x_between_l_and_h(T x) const { return (std::min(xl(), xh()) <= x && x <= std::max(xl(), xh())); }

      void set_cache_L_to_cache_x() { cache_xl_ = cache_x_; }
      void set_cache_H_to_cache_x() { cache_xh_ = cache_x_; }

      T x() const { return cache_x_.x(); }
      T f0() const { return cache_x_.f0(); }
      T xl() const { return cache_xl_.x(); }
      T f0_l() const { return cache_xl_.f0(); }
      T f1_l() const { return cache_xl_.f1(); }
      T xh() const { return cache_xh_.x(); }
      T f0_h() const { return cache_xh_.f0(); }
      T f1_h() const { return cache_xh_.f1(); }
      T dx_pp() const { return dx_history_.dx_pp(); }
      std::uintmax_t count() const { return count_; }

      F f_;
      detail::CacheOrder1<T> cache_x_;
      detail::CacheOrder1<T> cache_xl_;
      detail::CacheOrder1<T> cache_xh_;
      std::uintmax_t count_;
      T factor_;
      StepSizeHistory dx_history_;
   };
}  // namespace detail

template <class F, class T>
T newton_raphson_iterate(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::cout << "================= newton_raphson_iterate ====================" << std::endl;
   return detail::BracketHelper<F, T, detail::StepNewton<T>>::solve(f, x, xl, xh, digits, max_iter);
}

template <class F, class T>
inline T newton_raphson_iterate(F f, T guess, T min, T max, int digits) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   std::uintmax_t m = (std::numeric_limits<std::uintmax_t>::max)();
   return newton_raphson_iterate(f, guess, min, max, digits, m);
}

#endif 



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
