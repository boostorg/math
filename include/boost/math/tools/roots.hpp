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

#include <array>
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


namespace detail {
   ////// 1D Newton root finder explanation. //////
   //
   // Consider the following three root finding problems:
   //
   //   Case 1:  Find a root of the function f(x) given bounds [x_l_orig, x_h_orig] and
   //            an initial evaluation (x_0, f(x_0), f'(x_0)).
   //   Case 2U: Find a root of the function f(x) given evaluations at a low point
   //            (x_l, f(x_l), f'(x_l)) and a high point (x_h, f(x_h), f'(x_h)), where
   //            and f(x_l) and f(x_h) have the same sign, but have slopes that indicate
   //            that the function is sloping towards zero between x_l and x_h, forming
   //            a U shape (hence the name). The U may or may not cross zero.
   //   Case 2B: Find a root of the function f(x) given evaluations at a low point
   //            (x_l, f(x_l), f'(x_l)) and a high point (x_h, f(x_h), f'(x_h)), where
   //            and f(x_l) and f(x_h) have opposite signs that BRACKET a root.
   //
   // All three cases operate on similar data, but require different algorithms to
   // solve. These cases can also transitions between each other. Specifically...
   //
   //   Case 1:  The problem will remain in Case 1 as long as dx does not change sign.
   //            If dx changes sign, the problem will transition to Case 2U or Case 2B.
   //   Case 2U: The problem will transition to Case 2B if the criteria for Case 2B is
   //            met, i.e., f(x_l) and f(x_h) have opposite signs.
   //   Case 2B: Does not transition to other cases.
   //
   //
   ////// Code organization. //////
   //
   // The code for the root finding problem is divided into two parts: a state that
   // contains data and solvers that are dataless. State is stored in the Root1D_State
   // class. The following solver classes operate on the state:
   //      - SolverCase1
   //      - SolverCase2U
   //      - SolverCase2B
   // The process of finding a root begins by iterating with the solver for Case 1.
   // Iterations continue until one of four following things happen:
   //   - Returns success (f(x) is zero, ... etc)
   //   - Transition to Case 2U
   //   - Transition to Case 2B
   //   - Evaluated limit without being able to transition to Case 2U or Case 2B
   // Similarly, when iterating in Case 2U, iterations will continue until one of the
   // following happens:
   //   - Transition to Case 2B
   //   - Returns success (f(x) is zero, ... etc)
   //   - Returns failure (converging to non-zero extrema)
   // When iterating in Case 2B, iterations will continue until:
   //    - Returns success (f(x) is zero, ... etc)
   //
   // Enum for the possible cases for the root finding state machine
   enum class CaseFlag { case_1, case_2U, case_2B, success, failure };

   // Stores x, f(x), and f'(x) from a function evaluation
   template <class T, class Step>
   class Data_X01 {
   public:
      template <class F>
      Data_X01(F f, T x) : x_(x) { detail::unpack_tuple(f(x_), f0_, f1_); }

      // Store just a constant x value
      Data_X01(T x) : x_(x), f0_(static_cast<T>(1) / 0), f1_(static_cast<T>(1) / 0) {}

      T calc_dx() const { return Step().calc_dx(*this); }
      T calc_dx_newton() const { return -f0_ / f1_; }

      // This evaluation is a high bound if the Newton step is negative. This function checks the following
      // without division: -f(x) / f'(x) < 0
      bool is_bound_h() const { return sign(f0_) == sign(f1_); }

      // Returns true if *this and z BRACKET a root. For example between: f0(x_this) = 2 and f0(x_z) = -1 OR
      //                                                                  f0(x_this) = 0 and f0(x_z) = -1
      bool is_B(const Data_X01<T, Step>& z) const { return !is_same_sign_f0(z); }

      bool is_same_sign_f0(const Data_X01<T, Step>& z) const { return sign(f0_) == sign(z.f0()); }

      CaseFlag identify_status_flag(const Data_X01<T, Step>& h) const {
         if (is_B(h)) { return CaseFlag::case_2B; }
         if (is_U(h)) { return CaseFlag::case_2U; }
         return CaseFlag::failure;
      }

#ifdef BOOST_MATH_INSTRUMENT
      // Makes debugging easier if required
      void print() const { std::cout << "x: " << x_ << ", f0: " << f0_ << ", f1: " << f1_ << std::endl; }
#endif

      T x() const { return x_; }
      T f0() const { return f0_; }  // f(x)
      T f1() const { return f1_; }  // f'(x)

   private:
      bool is_sorted_x(const Data_X01<T, Step>& h) const { return x() <= h.x(); }

      // Returns true if *this and h form a U shape.
      bool is_U(const Data_X01<T, Step>& h) const {
         const auto& l = *this;
         if (!l.is_sorted_x(h)) { return h.is_U(l); }
         return (sign(l.f0()) * sign(l.f1()) == -1) &&  // Above zero and sloping down OR below zero and sloping up
                (sign(h.f0()) * sign(h.f1()) == +1);    // Above zero and sloping up OR below zero and sloping down
      }

      T x_;
      T f0_;  // f(x)
      T f1_;  // f'(x)
   };

   // Stores x and dx for a potential next root finding evaluation
   template <class T, class Step>
   class X_DX {
   public:
      X_DX(T x, T dx) : x_(x), dx_(dx) {}

      T x() const { return x_; }
      T dx() const { return dx_; }

   private:
      T x_;
      T dx_;
   };

   // Calculates the step with Newton's method
   template <class T>
   class StepNewton {
   public:
      T calc_dx(const Data_X01<T, StepNewton>& z) const { return -z.f0() / z.f1(); }
   };

   // The state of the 1D root finding problem. This state is acted upon by one of
   // several solvers
   template <class F, class T, class Step>
   class Root1D_State {
   private:
      enum class SideSelect { sign_f0,   // Updates based on the sign of f(x)
                              sign_dx,   // Updates based on the sign of dx
                              not_last,  // Updates entry not updated last
                              l,         // Updates the low bound
                              h };       // Updates the high bound

      // Allows template dispatch of update_eval
      template <SideSelect ValType>
      class Val {};

   public:
      Root1D_State(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter)
         : data_(xl, xh)
         , flag_(CaseFlag::case_1)
         , evaluator_(f, digits, max_iter)
      {
         update_eval_sign_dx(x);
      }

      // Syntactic sugar for type dispatch of update_eval
      template <typename... Args>
      void update_eval_sign_f0(Args... args) { update_eval(Val<SideSelect::sign_f0>(), args...); }
      template <typename... Args>
      void update_eval_sign_dx(Args... args) { update_eval(Val<SideSelect::sign_dx>(), args...); }
      template <typename... Args>
      void update_eval_not_last(Args... args) { update_eval(Val<SideSelect::not_last>(), args...); }
      template <typename... Args>
      void update_eval_l(Args... args) { update_eval(Val<SideSelect::l>(), args...); }
      template <typename... Args>
      void update_eval_h(Args... args) { update_eval(Val<SideSelect::h>(), args...); }

      // Filter arguments into preferred format
      template <SideSelect S>
      void update_eval(Val<S> vs, const T x) { update_eval(vs, calc_x_dx(x)); }
      template <SideSelect S>
      void update_eval(Val<S> vs, const X_DX<T, Step>& x_dx) { update_eval(vs, eval_next(x_dx), x_dx.dx()); }
      template <SideSelect S>
      void update_eval(Val<S> vs, const Data_X01<T, Step>& cx) { update_eval(vs, cx, calc_dx_true(cx)); }

      // Resolve SideSelect
      void update_eval(Val<SideSelect::sign_f0>, const Data_X01<T, Step>& cx, const T dx) {
         const bool is_bound_h = cx.is_same_sign_f0(eval_h());  // Updates high bound if f0 is same sign as f0_h
         update_eval(is_bound_h, cx, dx);
      }
      void update_eval(Val<SideSelect::sign_dx>, const Data_X01<T, Step>& cx, const T dx) {
         const bool is_bound_h = cx.is_bound_h();  // Updates high bound if dx is negative
         update_eval(is_bound_h, cx, dx);
      }
      void update_eval(Val<SideSelect::not_last>, const Data_X01<T, Step>& cx, const T dx) {
         const bool is_bound_h = !data_.ind_last();  // Updates high bound if last eval was not high
         update_eval(is_bound_h, cx, dx);
      }
      void update_eval(Val<SideSelect::l>, const Data_X01<T, Step>& cx, const T dx) {
         update_eval(false, cx, dx);  // Updates low bound
      }
      void update_eval(Val<SideSelect::h>, const Data_X01<T, Step>& cx, const T dx) {
         update_eval(true, cx, dx);  // Updates high bound
      }

      // Update the eval
      void update_eval(const bool is_bound_h, const Data_X01<T, Step>& cx, const T dx) {
         dx_history_.push(dx);
         data_.set(cx, is_bound_h);
         if (cx.f0() == 0) { flag_ = CaseFlag::success; }
      }

      void update_case_flag() {
         if (flag() == CaseFlag::success) { return; }  // Don't update if already successful

         // Find status flag if both evals are valid
         if (num_valid_evals() == 2) {
            flag_ = eval_l().identify_status_flag(eval_h());
         }
      }

      T midpoint() const { return (eval_l().x() + eval_h().x()) / 2; }
      bool is_inbounds(T x) const { return (x_l() <= x) && (x <= x_h()); }

      // NOTE: evals always contain valid x data, as both evals are initalized with the
      //       bounds of the root finding problem.
      T x_l() const { return eval_l().x(); }
      T x_h() const { return eval_h().x(); }

      // Returns the last eval that was evaluated
      const Data_X01<T, Step>& eval_last() const { return data_.eval_last(); }

      int num_x_only() const { return data_.num_x_only(); }
      int num_valid_evals() const { return 2 - num_x_only(); }

      // Extrapolates the last evaluation to calculate a next x value, and the distance dx
      // between the last x and this next x.
      X_DX<T, Step> extrapolate_last_to_calc_x_dx() const {
         const auto& cx = eval_last();
         const T x_next = cx.x() + cx.calc_dx();
         return calc_x_dx(x_next);
         // NOTE: the code below is incorrect because (x + dx) - x != dx
         // const T x = cx.x();
         // const T dx = cx.calc_dx();
         // return X_DX<T, Step>(x + dx, dx);
      }
      X_DX<T, Step> calc_x_dx(T x) const { return X_DX<T, Step>(x, calc_dx_true(x)); }

      // Calculates dx by evaluating x_next - x_prev. This is the actual floating point dx.
      T calc_dx_true(const Data_X01<T, Step>& cx) const { return calc_dx_true(cx.x()); }
      T calc_dx_true(const T& x) const { return x - eval_last().x(); }

#ifdef BOOST_MATH_INSTRUMENT
      // Makes debugging easier if required
      void print_status_flag() const {
         switch (flag()) {
            case CaseFlag::case_1: std::cout << "case_1" << std::endl; break;
            case CaseFlag::case_2U: std::cout << "case_2U" << std::endl; break;
            case CaseFlag::case_2B: std::cout << "case_2B" << std::endl; break;
            case CaseFlag::success: std::cout << "success" << std::endl; break;
            case CaseFlag::failure: std::cout << "failure" << std::endl; break;
         }
      }
#endif

      // Returns true if the two evals bracket a root where the sign of f(x) changes
      bool is_B() const { return eval_l().is_B(eval_h()); }

      void set_flag(CaseFlag flag) { flag_ = flag; }
      Data_X01<T, Step> eval_next(const X_DX<T, Step>& x_dx) { return evaluator_(x_dx); }
      void reset_dx_history() { dx_history_.reset(); }

      std::uintmax_t count() const { return evaluator_.count(); }
      T factor() const { return evaluator_.factor(); }
      T dx_pp() const { return dx_history_.dx_pp(); }
      std::uintmax_t num_fn_evals() const { return evaluator_.num_fn_evals(); }

      T get_active_bound() const { return data_.get_bound_from_unused_eval(); }

      const Data_X01<T, Step>& eval_l() const { return data_.eval_l(); }
      const Data_X01<T, Step>& eval_h() const { return data_.eval_h(); }
      CaseFlag flag() const { return flag_; }

   private:
      // Holds Evaluation Data
      class Root1D_Data {
      public:
         // Creates two evals whose x values are initialized with the bounds of the root finding
         // problem, but whose f0 and f1 values are invalid.
         Root1D_Data(T xl, T xh)
            : v_eval_({Data_X01<T, Step>(xl), Data_X01<T, Step>(xh)})
            , v_is_x_only_({{true, true}})
            , ind_last_(true) {}

         // Updates the eval using the bool is_h as an index.
         void set(const Data_X01<T, Step>& cx, bool is_h) {
            v_eval_[is_h] = cx;
            v_is_x_only_[is_h] = false;
            ind_last_ = is_h;
         }

         // Returns the x value of the only eval that contains an original bound
         T get_bound_from_unused_eval() const {
            BOOST_MATH_ASSERT(num_x_only() == 1);
            return v_is_x_only_[0] ? eval_l().x() : eval_h().x();
         }

         const Data_X01<T, Step>& eval_last() const { return v_eval_[ind_last_]; }
         const Data_X01<T, Step>& eval_l() const { return v_eval_[0]; }
         const Data_X01<T, Step>& eval_h() const { return v_eval_[1]; }
         
         int num_x_only() const { return v_is_x_only_[0] + v_is_x_only_[1]; }

         bool ind_last() const { return ind_last_; }  // Last evaluated index

      private:
         std::array<Data_X01<T, Step>, 2> v_eval_;
         std::array<bool, 2> v_is_x_only_;
         bool ind_last_;
      };

      // Stores the size of the last two steps. Calling reset will set all step sizes to infinity.
      class StepSizeHistory {
         public:
            StepSizeHistory() { reset(); }

            void reset() {
               // Need to static cast double infinity because for T = boost::math::concepts::real_concept, 
               // std::numeric_limits<T>::infinity() evaluates to 0.
               dx_p_ = dx_pp_ = static_cast<T>(std::numeric_limits<double>::infinity());
            }

            void push(T input) {
               dx_pp_ = dx_p_;
               dx_p_ = input;
            }

            T dx_pp() const { return dx_pp_; }

         private:
            T dx_p_;
            T dx_pp_;
      };

      // Stores the problem setup for the 1D root finding algorithm
      class Root1D_Evaluator {
      public:
         Root1D_Evaluator(F f, int digits, std::uintmax_t max_iter)
            : f_(f)
            , factor_(static_cast<T>(ldexp(1.0, 1 - digits)))
            , max_iter_(max_iter)
            , count_(max_iter) {}

         Data_X01<T, Step> operator()(const X_DX<T, Step>& x_dx) {
            if (count_ <= 0) {
               static const char* function = "boost::math::tools::detail::Root1D_Evaluator<%1%>";
               policies::raise_evaluation_error(function, "Ran out of iterations when finding root in boost::math::tools::detail::Root1D_Evaluator, attempted to evalutate %1%", x_dx.x(), boost::math::policies::policy<>());
            }
            --count_;
            return Data_X01<T, Step>(f_, x_dx.x());
         }

         std::uintmax_t num_fn_evals() const { return max_iter_ - count_; }

         T factor() const { return factor_; }
         std::uintmax_t count() const { return count_; }

      private:
         F f_;
         T factor_;
         std::uintmax_t max_iter_;
         std::uintmax_t count_;
      };

      // // // Data members // // //
      Root1D_Data data_;
      CaseFlag flag_;
      StepSizeHistory dx_history_;
      Root1D_Evaluator evaluator_;
   };

   // The abstract base class for the 1D root finding algorithm
   template <class F, class T, class Step>
   class SolverBase {
   public:
      // Solves the root finding problem for this status flag
      T solve(Root1D_State<F, T, Step>& b) const {
         X_DX<T, Step> x_dx = b.extrapolate_last_to_calc_x_dx();

         while (is_solver_current(b)) {
            // Calculate next x and dx
            x_dx = b.extrapolate_last_to_calc_x_dx();
            if (is_need_bisection(b, x_dx)) {
               b.reset_dx_history();
               x_dx = calc_next_bisection(b);
            }

            // Maybe exit
            if (is_exit_case_now(b, x_dx)) { break; }

            // Update eval
            this->update_eval(b, x_dx);
         }

         return x_dx.x();
      }

      bool is_need_bisection(Root1D_State<F, T, Step>& b, const X_DX<T, Step>& x_dx) const {
         using std::fabs;
         return !(b.is_inbounds(x_dx.x())) ||
                  fabs(b.dx_pp()) < fabs(x_dx.dx() * 4);  // Converges too slow
      }

      bool is_dx_small(const Root1D_State<F, T, Step>& b, const X_DX<T, Step> x_dx) const { 
         using std::fabs;
         return fabs(x_dx.dx()) <= fabs(x_dx.x()) * b.factor();
      }

      bool is_solver_current(const Root1D_State<F, T, Step>& b) const { return b.flag() == flag_this(); }

      virtual bool is_exit_case_now(Root1D_State<F, T, Step>& b, const X_DX<T, Step> x_dx) const = 0;
      virtual X_DX<T, Step> calc_next_bisection(Root1D_State<F, T, Step>& b) const = 0;
      virtual void update_eval(Root1D_State<F, T, Step>& b, const X_DX<T, Step>& x_dx) const = 0;
      virtual CaseFlag flag_this() const = 0;
   };

   // The solver for the case where only a single evulation (either low or high) exists.
   template <class F, class T, class Step>
   class SolverCase1 : public SolverBase<F, T, Step> {
   public:
      X_DX<T, Step> calc_next_bisection(Root1D_State<F, T, Step>& b) const override {
         const T x_limit = b.get_active_bound();
         const auto next_limit = b.calc_x_dx(x_limit);
         b.update_eval_not_last(next_limit);
         b.update_case_flag();
         return next_limit;
      }

      void update_eval(Root1D_State<F, T, Step>&b, const X_DX<T, Step>& x_dx) const override {
         b.update_eval_sign_dx(x_dx);
         b.update_case_flag();
      }

      bool is_exit_case_now(Root1D_State<F, T, Step>& b, const X_DX<T, Step> x_dx) const override {
         if (this->is_dx_small(b, x_dx)) {
            b.set_flag(CaseFlag::success); return true; }
         return b.flag() != flag_this();
      }

      CaseFlag flag_this() const override { return CaseFlag::case_1; }
   };

   // Base class for SolverCase2U and SolverCase2B
   template <class F, class T, class Step>
   class SolverCase2Base : public SolverBase<F, T, Step> {
   public:
      X_DX<T, Step> calc_next_bisection(Root1D_State<F, T, Step>& b) const override {
         const auto x_mid = b.midpoint();
         if (x_mid == b.x_l() || x_mid == b.x_h()) {
            b.set_flag(this->flag_stall());
         }
         return b.calc_x_dx(x_mid);
      }

      virtual CaseFlag flag_stall() const = 0;
   };

   // The solver for Case 2U
   template <class F, class T, class Step>
   class SolverCase2U : public SolverCase2Base<F, T, Step> {
   public:
      void update_eval(Root1D_State<F, T, Step>& b, const X_DX<T, Step>& x_dx) const override {
         const auto cx = b.eval_next(x_dx);

         // Transition to Case 2B because will bracket a root with the new evaluation
         if (b.eval_last().is_B(cx)) {
            b.set_flag(CaseFlag::case_2B);

         // If the newton step is increasing this suggests that the solver is converging toward a
         // local minimum that is not a root. This suggests failure.
         } else if (is_newton_dx_bigger_than_prev(b, cx)) {
            b.set_flag(CaseFlag::failure);
         }

         b.update_eval_sign_dx(cx, x_dx.dx());
      }

      bool is_newton_dx_bigger_than_prev(Root1D_State<F, T, Step>& b, const Data_X01<T, Step>& cx) const {
         using std::fabs;
         // NOTE: b.eval_last() isn't generally on the same side of the root as cx
         const auto& cx_prev_side_eval = cx.is_bound_h() ? b.eval_h() : b.eval_l();
         return fabs(cx_prev_side_eval.calc_dx_newton()) < fabs(cx.calc_dx_newton());  // Uses Newton step
      }

      bool is_exit_case_now(Root1D_State<F, T, Step>& b, const X_DX<T, Step> x_dx) const override {
         if (x_dx.dx() == 0) {
            b.set_flag(CaseFlag::failure); return true;
         }
         return b.flag() != flag_this();
      }

      CaseFlag flag_this() const override { return CaseFlag::case_2U; }
      CaseFlag flag_stall() const override { return CaseFlag::failure; }
   };

   // The solver for Case 2B
   template <class F, class T, class Step>
   class SolverCase2B : public SolverCase2Base<F, T, Step> {
   public:
      void update_eval(Root1D_State<F, T, Step>& b, const X_DX<T, Step>& x_dx) const override {
         b.update_eval_sign_f0(x_dx);
      }

      bool is_exit_case_now(Root1D_State<F, T, Step>& b, const X_DX<T, Step> x_dx) const override {
         if (this->is_dx_small(b, x_dx)) {
            b.set_flag(CaseFlag::success); return true; }
         return b.flag() != flag_this();
      }

      CaseFlag flag_this() const override { return CaseFlag::case_2B; }
      CaseFlag flag_stall() const override { return CaseFlag::success; }
   };

   // Iterates on the root finding problem using the solver associated with the current case flag.
   // Returns the result.
   template <class F, class T, class Step>
   std::pair<bool, T> solve_for_state(Root1D_State<F, T, Step>& state, T x_final) {
      switch (state.flag()) {
         case CaseFlag::case_1: return solve_for_state(state, SolverCase1<F, T, Step>().solve(state));
         case CaseFlag::case_2B: return solve_for_state(state, SolverCase2B<F, T, Step>().solve(state));
         case CaseFlag::case_2U: return solve_for_state(state, SolverCase2U<F, T, Step>().solve(state));
         case CaseFlag::failure: return std::make_pair(false, 0);
         case CaseFlag::success: return std::make_pair(true, x_final);
      }
      return std::make_pair(false, 0);  // Unreachable but suppresses compiler warning
   }
   template <class F, class T, class Step>
   std::pair<bool, T> solve_for_state(Root1D_State<F, T, Step>& state) {
      return solve_for_state(state, state.eval_last().x());
   }

   // Makes three attempts to find a root:
   //   1. Local minimization starting with the initial x
   //   2. Bracketing the root in the direction of the initial dx if there is a sign flip
   //      between f(x) and f(x_bound in direction of dx)
   //   3. Bracketing the root in the opposite direction of the initial dx if there is a
   //      sign flip between f(x) and f(x_bound in opposite direction of dx).
   //
   // If all three attempts fail, an error is thrown.
   //
   //
   template <class F, class T, class Step>
   T solve_3_attempts(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>()))) {
      // Create Root1D_State
      Root1D_State<F, T, Step> state(f, x, xl, xh, digits, max_iter);

      // Get CacheX
      const auto eval_x_orig = state.eval_last();

      auto fn_form_bracket_with_side = [&](const bool ind_bound){
         if (ind_bound) {
            // NOTE: These two lines cannot be swapped in case xh is a solution. This is tested.
            state.update_eval_l(eval_x_orig);  // First to original eval
            state.update_eval_h(xh);           // Then set bound
         } else {
            // NOTE: These two lines cannot be swapped in case xl is a solution. This is tested.
            state.update_eval_h(eval_x_orig);  // First to original eval
            state.update_eval_l(xl);           // Then set bound
         }

         // state.reset_dx_history();  // Not required for correctness
         state.update_case_flag();  // Sets to case_2U, case_2B, or failure
      };

      auto fn_return = [&](const Root1D_State<F, T, Step>& s, const std::pair<bool, T> p) {
         max_iter = s.num_fn_evals();  // Update max_iter
         return p.second;
      };

      // Solve problem
      const auto p = solve_for_state(state);

      // If local minimization was successful, return root
      if (p.first) { return fn_return(state, p); }

      // Attempt to bracket root in the direction of the initial guess
      fn_form_bracket_with_side(!eval_x_orig.is_bound_h());
      if (state.is_B()) { return fn_return(state, solve_for_state(state)); }

      // Attempt to bracket root in the opposite direction of the initial guess
      fn_form_bracket_with_side(eval_x_orig.is_bound_h());
      if (state.is_B()) { return fn_return(state, solve_for_state(state)); }

      static const char* function = "boost::math::tools::detail::solve<%3%>";
      const T x_last = state.eval_last().x();
      return policies::raise_evaluation_error(function, "Unable to bracket root in boost::math::tools::detail::solve, last closest guess was %1%", x_last, boost::math::policies::policy<>());
   }
}  // namespace detail

template <class F, class T>
T newton_raphson_iterate(F f, T x, T xl, T xh, int digits, std::uintmax_t& max_iter) noexcept(policies::is_noexcept_error_policy<policies::policy<> >::value&& BOOST_MATH_IS_FLOAT(T) && noexcept(std::declval<F>()(std::declval<T>())))
{
   return detail::solve_3_attempts<F, T, detail::StepNewton<T>>(f, x, xl, xh, digits, max_iter);
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
