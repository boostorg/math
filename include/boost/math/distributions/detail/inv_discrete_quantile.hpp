//  Copyright John Maddock 2007.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_DETAIL_INV_DISCRETE_QUANTILE
#define BOOST_MATH_DISTRIBUTIONS_DETAIL_INV_DISCRETE_QUANTILE

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/cstdint.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#include <boost/math/tools/tuple.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace detail{

//
// Functor for root finding algorithm:
//
template <class Dist>
struct distribution_quantile_finder
{
   typedef typename Dist::value_type value_type;
   typedef typename Dist::policy_type policy_type;

   BOOST_MATH_GPU_ENABLED distribution_quantile_finder(const Dist d, value_type p, bool c)
      : dist(d), target(p), comp(c) {}

   BOOST_MATH_GPU_ENABLED value_type operator()(value_type const& x)
   {
      return comp ? value_type(target - cdf(complement(dist, x))) : value_type(cdf(dist, x) - target);
   }

private:
   Dist dist;
   value_type target;
   bool comp;
};
//
// The purpose of adjust_bounds, is to toggle the last bit of the
// range so that both ends round to the same integer, if possible.
// If they do both round the same then we terminate the search
// for the root *very* quickly when finding an integer result.
// At the point that this function is called we know that "a" is
// below the root and "b" above it, so this change can not result
// in the root no longer being bracketed.
//
template <class Real, class Tol>
BOOST_MATH_GPU_ENABLED void adjust_bounds(Real& /* a */, Real& /* b */, Tol const& /* tol */){}

template <class Real>
BOOST_MATH_GPU_ENABLED void adjust_bounds(Real& /* a */, Real& b, tools::equal_floor const& /* tol */)
{
   BOOST_MATH_STD_USING
   b -= tools::epsilon<Real>() * b;
}

template <class Real>
BOOST_MATH_GPU_ENABLED void adjust_bounds(Real& a, Real& /* b */, tools::equal_ceil const& /* tol */)
{
   BOOST_MATH_STD_USING
   a += tools::epsilon<Real>() * a;
}

template <class Real>
BOOST_MATH_GPU_ENABLED void adjust_bounds(Real& a, Real& b, tools::equal_nearest_integer const& /* tol */)
{
   BOOST_MATH_STD_USING
   a += tools::epsilon<Real>() * a;
   b -= tools::epsilon<Real>() * b;
}
//
// This is where all the work is done:
//
template <class Dist, class Tolerance>
BOOST_MATH_GPU_ENABLED typename Dist::value_type 
   do_inverse_discrete_quantile(
      const Dist& dist,
      const typename Dist::value_type& p,
      bool comp,
      typename Dist::value_type guess,
      const typename Dist::value_type& multiplier,
      typename Dist::value_type adder,
      const Tolerance& tol,
      BOOST_MATH_NAMESPACE::uintmax_t& max_iter)
{
   typedef typename Dist::value_type value_type;
   typedef typename Dist::policy_type policy_type;

   constexpr auto function = "boost::math::do_inverse_discrete_quantile<%1%>";

   BOOST_MATH_STD_USING

   distribution_quantile_finder<Dist> f(dist, p, comp);
   //
   // Max bounds of the distribution:
   //
   value_type min_bound, max_bound;
   BOOST_MATH_NAMESPACE::tie(min_bound, max_bound) = support(dist);

   if(guess > max_bound)
      guess = max_bound;
   if(guess < min_bound)
      guess = min_bound;

   value_type fa = f(guess);
   BOOST_MATH_NAMESPACE::uintmax_t count = max_iter - 1;
   value_type fb(fa), a(guess), b =0; // Compiler warning C4701: potentially uninitialized local variable 'b' used

   if(fa == 0)
      return guess;

   //
   // For small expected results, just use a linear search:
   //
   if(guess < 10)
   {
      b = a;
      while((a < 10) && (fa * fb >= 0))
      {
         if(fb <= 0)
         {
            a = b;
            b = a + 1;
            if(b > max_bound)
               b = max_bound;
            fb = f(b);
            --count;
            if(fb == 0)
               return b;
            if(a == b)
               return b; // can't go any higher!
         }
         else
         {
            b = a;
            a = BOOST_MATH_GPU_SAFE_MAX(value_type(b - 1), value_type(0));
            if(a < min_bound)
               a = min_bound;
            fa = f(a);
            --count;
            if(fa == 0)
               return a;
            if(a == b)
               return a;  //  We can't go any lower than this!
         }
      }
   }
   //
   // Try and bracket using a couple of additions first, 
   // we're assuming that "guess" is likely to be accurate
   // to the nearest int or so:
   //
   else if((adder != 0) && (a + adder != a))
   {
      //
      // If we're looking for a large result, then bump "adder" up
      // by a bit to increase our chances of bracketing the root:
      //
      //adder = BOOST_MATH_GPU_SAFE_MAX(adder, 0.001f * guess);
      if(fa < 0)
      {
         b = a + adder;
         if(b > max_bound)
            b = max_bound;
      }
      else
      {
         b = BOOST_MATH_GPU_SAFE_MAX(value_type(a - adder), value_type(0));
         if(b < min_bound)
            b = min_bound;
      }
      fb = f(b);
      --count;
      if(fb == 0)
         return b;
      if(count && (fa * fb >= 0))
      {
         //
         // We didn't bracket the root, try 
         // once more:
         //
         a = b;
         fa = fb;
         if(fa < 0)
         {
            b = a + adder;
            if(b > max_bound)
               b = max_bound;
         }
         else
         {
            b = BOOST_MATH_GPU_SAFE_MAX(value_type(a - adder), value_type(0));
            if(b < min_bound)
               b = min_bound;
         }
         fb = f(b);
         --count;
      }
      if(a > b)
      {
         BOOST_MATH_GPU_SAFE_SWAP(a, b);
         BOOST_MATH_GPU_SAFE_SWAP(fa, fb);
      }
   }
   //
   // If the root hasn't been bracketed yet, try again
   // using the multiplier this time:
   //
   if((BOOST_MATH_NAMESPACE::sign)(fb) == (BOOST_MATH_NAMESPACE::sign)(fa))
   {
      if(fa < 0)
      {
         //
         // Zero is to the right of x2, so walk upwards
         // until we find it:
         //
         while(((BOOST_MATH_NAMESPACE::sign)(fb) == (BOOST_MATH_NAMESPACE::sign)(fa)) && (a != b))
         {
            if(count == 0)
               return policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", b, policy_type()); // LCOV_EXCL_LINE
            a = b;
            fa = fb;
            b *= multiplier;
            if(b > max_bound)
               b = max_bound;
            fb = f(b);
            --count;
            BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
         }
      }
      else
      {
         //
         // Zero is to the left of a, so walk downwards
         // until we find it:
         //
         while(((BOOST_MATH_NAMESPACE::sign)(fb) == (BOOST_MATH_NAMESPACE::sign)(fa)) && (a != b))
         {
            if(fabs(a) < tools::min_value<value_type>())
            {
               // Escape route just in case the answer is zero!
               max_iter -= count;
               max_iter += 1;
               return 0;
            }
            if(count == 0)
               return policies::raise_evaluation_error(function, "Unable to bracket root, last nearest value was %1%", a, policy_type()); // LCOV_EXCL_LINE
            b = a;
            fb = fa;
            a /= multiplier;
            if(a < min_bound)
               a = min_bound;
            fa = f(a);
            --count;
            BOOST_MATH_INSTRUMENT_CODE("a = " << a << " b = " << b << " fa = " << fa << " fb = " << fb << " count = " << count);
         }
      }
   }
   max_iter -= count;
   if(fa == 0)
      return a;
   if(fb == 0)
      return b;
   if(a == b)
      return b;  // Ran out of bounds trying to bracket - there is no answer!
   //
   // Adjust bounds so that if we're looking for an integer
   // result, then both ends round the same way:
   //
   adjust_bounds(a, b, tol);
   //
   // We don't want zero or denorm lower bounds:
   //
   if(a < tools::min_value<value_type>())
      a = tools::min_value<value_type>();
   //
   // Go ahead and find the root:
   //
   BOOST_MATH_NAMESPACE::pair<value_type, value_type> r = toms748_solve(f, a, b, fa, fb, tol, count, policy_type());
   max_iter += count;
   if (max_iter >= policies::get_max_root_iterations<policy_type>())
   {
      return policies::raise_evaluation_error<value_type>(function, "Unable to locate solution in a reasonable time:" // LCOV_EXCL_LINE
         " either there is no answer to quantile or the answer is infinite.  Current best guess is %1%", r.first, policy_type()); // LCOV_EXCL_LINE
   }
   BOOST_MATH_INSTRUMENT_CODE("max_iter = " << max_iter << " count = " << count);
   return (r.first + r.second) / 2;
}
//
// Rounding the real-valued root to an integer.
//
// The root finder above solves the *continuous* cdf(x) = p, and near an
// integer k the root it finds is only known to within the error of the cdf
// evaluated slightly off that integer.  Simply taking floor or ceil of it
// then gets the integer wrong by one whenever p is within a few ulps of
// cdf(k): https://github.com/boostorg/math/issues/935
//
// So we only use the real root as a starting point and settle the answer
// using the cdf at integers, which is what callers compare against:
//
//   round up:   the smallest k with cdf(k) >= p,
//   round down: the largest k with cdf(k) <= p.
//
// Round up is the usual definition of a discrete quantile (R, SciPy, Julia
// and Mathematica all use it).  Where a run of consecutive integers all have
// cdf(k) == p exactly (common as p -> 1) both take these definitions
// literally: round up gives the first of the run and round down the last.
// For the complement (c == true) read "cdf(k) - p" as
// "p - cdf(complement(d, k))", which also increases with k, so round up
// gives the smallest k with cdf(complement(d, k)) <= p.
//
// The support need not end on an integer (the binomial accepts fractional
// trials), so integers above it are handled without calling the cdf: there
// cdf(k) == 1 and cdf(complement(d, k)) == 0.
//
template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type discrete_quantile_residual(const Dist& d, const typename Dist::value_type& k, const typename Dist::value_type& p, bool c)
{
   typedef typename Dist::value_type value_type;
   if (k > support(d).second)
      return c ? p : value_type(1 - p);
   return c ? value_type(p - cdf(complement(d, k))) : value_type(cdf(d, k) - p);
}

//
// The smallest integer k in [lo, hi] at which the residual is >= 0 (or > 0
// when "strict"), searching from k.  The residual increases with k, and is
// taken to pass at hi without evaluating it there.  The answer is normally
// within one of k, but where the cdf rounds to the same value over a long
// run of integers (as p -> 1) it may be far away, so gallop, then bisect.
//
template <class Dist>
BOOST_MATH_GPU_ENABLED inline bool discrete_quantile_passes(const Dist& d, const typename Dist::value_type& k, const typename Dist::value_type& p, bool c, bool strict)
{
   typename Dist::value_type r = discrete_quantile_residual(d, k, p, c);
   return strict ? (r > 0) : (r >= 0);
}

template <class Dist>
BOOST_MATH_GPU_ENABLED typename Dist::value_type first_passing_integer(const Dist& d, typename Dist::value_type k, const typename Dist::value_type& p, bool c, bool strict, const typename Dist::value_type& lo, const typename Dist::value_type& hi)
{
   BOOST_MATH_STD_USING
   typedef typename Dist::value_type value_type;
   constexpr auto function = "boost::math::quantile<%1%>(%1%, %1%)";
   if (k < lo)
      k = lo;
   if (k > hi)
      k = hi;
   // Beyond 1/epsilon every value is an integer and k + 1 == k, so there is
   // nothing to settle:
   if (k + 1 == k)
      return k;
   // Each evaluation after the first counts against the root finder's
   // iteration budget, and running out of it is an error:
   BOOST_MATH_NAMESPACE::uintmax_t steps = policies::get_max_root_iterations<typename Dist::policy_type>();
   value_type pass = k, fail = k;  // Once both are known, fail < pass.
   value_type step = 1;
   if ((k == hi) || discrete_quantile_passes(d, k, p, c, strict))
   {
      for (;;)
      {
         if (pass == lo)
            return lo;
         if (steps-- == 0)
            return policies::raise_evaluation_error<value_type>(function, "Unable to settle the integer result, last value tried was %1%", pass, typename Dist::policy_type()); // LCOV_EXCL_LINE
         value_type t = (pass - lo > step) ? value_type(pass - step) : lo;
         if (!discrete_quantile_passes(d, t, p, c, strict))
         {
            fail = t;
            break;
         }
         pass = t;
         step *= 2;
      }
   }
   else
   {
      for (;;)
      {
         value_type t = (hi - fail > step) ? value_type(fail + step) : hi;
         if (t == hi)
         {
            pass = t;
            break;
         }
         if (steps-- == 0)
            return policies::raise_evaluation_error<value_type>(function, "Unable to settle the integer result, last value tried was %1%", fail, typename Dist::policy_type()); // LCOV_EXCL_LINE
         if (discrete_quantile_passes(d, t, p, c, strict))
         {
            pass = t;
            break;
         }
         fail = t;
         step *= 2;
      }
   }
   while (pass - fail > 1)
   {
      if (steps-- == 0)
         return policies::raise_evaluation_error<value_type>(function, "Unable to settle the integer result, last value tried was %1%", pass, typename Dist::policy_type()); // LCOV_EXCL_LINE
      value_type mid = fail + floor((pass - fail) / 2);
      if (discrete_quantile_passes(d, mid, p, c, strict))
         pass = mid;
      else
         fail = mid;
   }
   return pass;
}

//
// Integers only: the last one in the support is floor(support.second), and
// the one after it (where cdf == 1, so the residual is positive unless
// p == 1) bounds the searches.
//
template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_to_ceil(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c)
{
   BOOST_MATH_STD_USING
   typedef typename Dist::value_type value_type;
   const value_type lo = ceil(support(d).first);
   const value_type hi = floor(support(d).second) + 1;
   return first_passing_integer(d, value_type(ceil(result)), p, c, false, lo, hi);
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_to_floor(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c)
{
   BOOST_MATH_STD_USING
   typedef typename Dist::value_type value_type;
   const value_type lo = ceil(support(d).first);
   const value_type hi = floor(support(d).second) + 1;
   // One below the first k with cdf(k) > p.  If even cdf(lo) > p there is no
   // answer, and lo is the best we can do:
   const value_type k = first_passing_integer(d, value_type(floor(result) + 1), p, c, true, lo, hi);
   return k > lo ? value_type(k - 1) : lo;
}

//
// integer_round_nearest rounds result + 0.5 down, but "nearest" is not a
// statement about the cdf at integers, so it keeps the original logic:
// prefer ceil(result) if it is an exact root, else floor(result), then move
// down to the smallest k with cdf(k) >= p (what round up gives) if that is
// lower.  So the answer is the nearest integer, but never above round up's.
//
template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_to_nearest_floor(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c)
{
   BOOST_MATH_STD_USING
   typedef typename Dist::value_type value_type;
   const value_type lo = ceil(support(d).first);
   const value_type hi = floor(support(d).second) + 1;
   value_type k = floor(result);
   if (k < lo)
      k = lo;
   if (k > hi)
      k = hi;
   const value_type cc = ceil(result);
   if ((cc != k) && (cc >= lo) && (cc <= hi) && (discrete_quantile_residual(d, cc, p, c) == 0))
      return first_passing_integer(d, cc, p, c, false, lo, cc);
   if (discrete_quantile_passes(d, k, p, c, false))
      return first_passing_integer(d, k, p, c, false, lo, k);
   return k;
}

//
// Rounding for distributions that find their quantile some other way (in
// closed form, or by summing the pdf) and only need the result settled on
// the integer the policy asks for.  "result" is the real-valued quantile,
// or an integer within a step or so of the answer.
//
template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_discrete_quantile(const Dist&, typename Dist::value_type result, typename Dist::value_type, bool, const policies::discrete_quantile<policies::real>&)
{
   return result;
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_discrete_quantile(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_up>&)
{
   return round_to_ceil(d, result, p, c);
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_discrete_quantile(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_down>&)
{
   return round_to_floor(d, result, p, c);
}

//
// integer_round_outwards rounds lower quantiles (probability below 0.5) down
// and upper ones up; integer_round_inwards the reverse.  Both are a choice of
// one of the other two policies:
//
template <class T>
BOOST_MATH_GPU_ENABLED inline bool is_lower_quantile(const T& p, bool c)
{
   return (c ? T(1 - p) : p) < 0.5f;
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_discrete_quantile(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_outwards>&)
{
   return is_lower_quantile(p, c)
      ? round_discrete_quantile(d, result, p, c, policies::discrete_quantile<policies::integer_round_down>())
      : round_discrete_quantile(d, result, p, c, policies::discrete_quantile<policies::integer_round_up>());
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_discrete_quantile(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_inwards>&)
{
   return is_lower_quantile(p, c)
      ? round_discrete_quantile(d, result, p, c, policies::discrete_quantile<policies::integer_round_up>())
      : round_discrete_quantile(d, result, p, c, policies::discrete_quantile<policies::integer_round_down>());
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_discrete_quantile(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c, const policies::discrete_quantile<policies::integer_round_nearest>&)
{
   return round_to_nearest_floor(d, typename Dist::value_type(result + 0.5f), p, c);
}
//
// Now finally are the public API functions.
// There is one overload for each policy,
// each one is responsible for selecting the correct
// termination condition, and rounding the result
// to an int where required.
//
template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type 
   inverse_discrete_quantile(
      const Dist& dist,
      typename Dist::value_type p,
      bool c,
      const typename Dist::value_type& guess,
      const typename Dist::value_type& multiplier,
      const typename Dist::value_type& adder,
      const policies::discrete_quantile<policies::real>&,
      BOOST_MATH_NAMESPACE::uintmax_t& max_iter)
{
   if(p > 0.5)
   {
      p = 1 - p;
      c = !c;
   }
   if(discrete_quantile_residual(dist, typename Dist::value_type(0), p, c) >= 0)
      return 0;  // cdf(0) >= p already: see round_to_ceil / round_to_floor.
   return do_inverse_discrete_quantile(
      dist, 
      p, 
      c,
      guess, 
      multiplier, 
      adder, 
      tools::eps_tolerance<typename Dist::value_type>(policies::digits<typename Dist::value_type, typename Dist::policy_type>()),
      max_iter);
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type 
   inverse_discrete_quantile(
      const Dist& dist,
      const typename Dist::value_type& p,
      bool c,
      const typename Dist::value_type& guess,
      const typename Dist::value_type& multiplier,
      const typename Dist::value_type& adder,
      const policies::discrete_quantile<policies::integer_round_down>&,
      BOOST_MATH_NAMESPACE::uintmax_t& max_iter)
{
   typedef typename Dist::value_type value_type;
   BOOST_MATH_STD_USING
   if(discrete_quantile_residual(dist, typename Dist::value_type(0), p, c) >= 0)
      return 0;  // cdf(0) >= p already: see round_to_ceil / round_to_floor.
   return round_to_floor(dist, do_inverse_discrete_quantile(
      dist, 
      p, 
      c,
      (guess < 1 ? value_type(1) : floor(guess)), 
      multiplier, 
      adder, 
      tools::equal_floor(),
      max_iter), p, c);
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type 
   inverse_discrete_quantile(
      const Dist& dist,
      const typename Dist::value_type& p,
      bool c,
      const typename Dist::value_type& guess,
      const typename Dist::value_type& multiplier,
      const typename Dist::value_type& adder,
      const policies::discrete_quantile<policies::integer_round_up>&,
      BOOST_MATH_NAMESPACE::uintmax_t& max_iter)
{
   BOOST_MATH_STD_USING
   if(discrete_quantile_residual(dist, typename Dist::value_type(0), p, c) >= 0)
      return 0;  // cdf(0) >= p already: see round_to_ceil / round_to_floor.
   return round_to_ceil(dist, do_inverse_discrete_quantile(
      dist, 
      p, 
      c,
      ceil(guess), 
      multiplier, 
      adder, 
      tools::equal_ceil(),
      max_iter), p, c);
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type 
   inverse_discrete_quantile(
      const Dist& dist,
      const typename Dist::value_type& p,
      bool c,
      const typename Dist::value_type& guess,
      const typename Dist::value_type& multiplier,
      const typename Dist::value_type& adder,
      const policies::discrete_quantile<policies::integer_round_nearest>&,
      BOOST_MATH_NAMESPACE::uintmax_t& max_iter)
{
   typedef typename Dist::value_type value_type;
   BOOST_MATH_STD_USING
   if(discrete_quantile_residual(dist, typename Dist::value_type(0), p, c) >= 0)
      return 0;  // cdf(0) >= p already: see round_to_ceil / round_to_floor.
   //
   // Note that we adjust the guess to the nearest half-integer:
   // this increase the chances that we will bracket the root
   // with two results that both round to the same integer quickly.
   //
   return round_to_nearest_floor(dist, do_inverse_discrete_quantile(
      dist, 
      p, 
      c,
      (guess < 0.5f ? value_type(1.5f) : floor(guess + 0.5f) + 0.5f), 
      multiplier, 
      adder, 
      tools::equal_nearest_integer(),
      max_iter) + 0.5f, p, c);
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type 
   inverse_discrete_quantile(
      const Dist& dist,
      const typename Dist::value_type& p,
      bool c,
      const typename Dist::value_type& guess,
      const typename Dist::value_type& multiplier,
      const typename Dist::value_type& adder,
      const policies::discrete_quantile<policies::integer_round_outwards>&,
      BOOST_MATH_NAMESPACE::uintmax_t& max_iter)
{
   // The up and down overloads also return 0 when cdf(0) >= p already.
   return is_lower_quantile(p, c)
      ? inverse_discrete_quantile(dist, p, c, guess, multiplier, adder, policies::discrete_quantile<policies::integer_round_down>(), max_iter)
      : inverse_discrete_quantile(dist, p, c, guess, multiplier, adder, policies::discrete_quantile<policies::integer_round_up>(), max_iter);
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type 
   inverse_discrete_quantile(
      const Dist& dist,
      const typename Dist::value_type& p,
      bool c,
      const typename Dist::value_type& guess,
      const typename Dist::value_type& multiplier,
      const typename Dist::value_type& adder,
      const policies::discrete_quantile<policies::integer_round_inwards>&,
      BOOST_MATH_NAMESPACE::uintmax_t& max_iter)
{
   // The up and down overloads also return 0 when cdf(0) >= p already.
   return is_lower_quantile(p, c)
      ? inverse_discrete_quantile(dist, p, c, guess, multiplier, adder, policies::discrete_quantile<policies::integer_round_up>(), max_iter)
      : inverse_discrete_quantile(dist, p, c, guess, multiplier, adder, policies::discrete_quantile<policies::integer_round_down>(), max_iter);
}

} BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_DISTRIBUTIONS_DETAIL_INV_DISCRETE_QUANTILE

