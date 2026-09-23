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

namespace boost{ namespace math{ namespace detail{

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
      boost::math::uintmax_t& max_iter)
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
   boost::math::tie(min_bound, max_bound) = support(dist);

   if(guess > max_bound)
      guess = max_bound;
   if(guess < min_bound)
      guess = min_bound;

   value_type fa = f(guess);
   boost::math::uintmax_t count = max_iter - 1;
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
   if((boost::math::sign)(fb) == (boost::math::sign)(fa))
   {
      if(fa < 0)
      {
         //
         // Zero is to the right of x2, so walk upwards
         // until we find it:
         //
         while(((boost::math::sign)(fb) == (boost::math::sign)(fa)) && (a != b))
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
         while(((boost::math::sign)(fb) == (boost::math::sign)(fa)) && (a != b))
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
   boost::math::pair<value_type, value_type> r = toms748_solve(f, a, b, fa, fb, tol, count, policy_type());
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
// When a run of consecutive integers all have cdf(k) == p exactly (common
// as p -> 1), round up takes the last of them and round down the first,
// as this code always has.  For the complement (c == true) read
// "cdf(k) - p" as "p - cdf(complement(d, k))", which also increases with k.
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

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_to_floor(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c)
{
   BOOST_MATH_STD_USING
   typedef typename Dist::value_type value_type;
   // Integers only: the last one in the support is floor(support.second).
   const value_type lo = ceil(support(d).first);
   const value_type hi = floor(support(d).second);
   value_type k = floor(result);
   if (k < lo)
      k = lo;
   if (k > hi)
      k = hi;
   value_type gk = discrete_quantile_residual(d, k, p, c);
   // Step down until cdf(k) <= p:
   while ((gk > 0) && (k > lo))
   {
      --k;
      gk = discrete_quantile_residual(d, k, p, c);
   }
   // Step up while the next integer still has cdf <= p.  If it hits a run
   // with cdf == p we want the first of that run, so stop there:
   while (k < hi)
   {
      value_type gn = discrete_quantile_residual(d, value_type(k + 1), p, c);
      if (gn > 0)
         break;
      if (gn == 0)
      {
         if (gk < 0)
            return k + 1; // First member of the run.
         break;           // Already inside the run.
      }
      ++k;
      gk = gn;
   }
   // If we started inside a run with cdf == p, move to its first member:
   if (gk == 0)
   {
      while (k > lo)
      {
         if (discrete_quantile_residual(d, value_type(k - 1), p, c) < 0)
            break;
         --k;
      }
   }
   return k;
}

template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_to_ceil(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c)
{
   BOOST_MATH_STD_USING
   typedef typename Dist::value_type value_type;
   // Integers only.  Rounding up may land one past the support (where
   // cdf == 1) when p exceeds the cdf at the last integer in it.
   const value_type lo = ceil(support(d).first);
   const value_type hi = floor(support(d).second) + 1;
   value_type k = ceil(result);
   if (k < lo)
      k = lo;
   if (k > hi)
      k = hi;
   value_type gk = discrete_quantile_residual(d, k, p, c);
   // Step up until cdf(k) >= p:
   while ((gk < 0) && (k < hi))
   {
      ++k;
      gk = discrete_quantile_residual(d, k, p, c);
   }
   // Step down while the previous integer still has cdf >= p.  If it hits
   // a run with cdf == p we want the last of that run, so stop there:
   while (k > lo)
   {
      value_type gp = discrete_quantile_residual(d, value_type(k - 1), p, c);
      if (gp < 0)
         break;
      if (gp == 0)
      {
         if (gk > 0)
            return k - 1; // Last member of the run.
         break;           // Already inside the run.
      }
      --k;
      gk = gp;
   }
   // If we started inside a run with cdf == p, move to its last member:
   if (gk == 0)
   {
      while (k < hi)
      {
         if (discrete_quantile_residual(d, value_type(k + 1), p, c) > 0)
            break;
         ++k;
      }
   }
   return k;
}

//
// integer_round_nearest rounds result + 0.5 down, but "nearest" is not a
// statement about the cdf at integers, so it keeps the original logic:
// prefer ceil(result) if it is an exact root, else floor(result), then
// move to the smallest integer that is still an exact root.
//
template <class Dist>
BOOST_MATH_GPU_ENABLED inline typename Dist::value_type round_to_nearest_floor(const Dist& d, typename Dist::value_type result, typename Dist::value_type p, bool c)
{
   BOOST_MATH_STD_USING
   typename Dist::value_type cc = ceil(result);
   typename Dist::value_type pp = cc <= support(d).second ? c ? cdf(complement(d, cc)) : cdf(d, cc) : 1;
   if(pp == p)
      result = cc;
   else
      result = floor(result);
   while(result != 0)
   {
      #ifdef BOOST_MATH_HAS_GPU_SUPPORT
      cc = floor(::nextafter(result, -tools::max_value<typename Dist::value_type>()));
      #else
      cc = floor(float_prior(result));
      #endif
      if(cc < support(d).first)
         break;
      pp = c ? cdf(complement(d, cc)) : cdf(d, cc);
      if(c ? pp > p : pp < p)
         break;
      result = cc;
   }

   return result;
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
      boost::math::uintmax_t& max_iter)
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
      const policies::discrete_quantile<policies::integer_round_outwards>&,
      boost::math::uintmax_t& max_iter)
{
   typedef typename Dist::value_type value_type;
   BOOST_MATH_STD_USING
   typename Dist::value_type pp = c ? 1 - p : p;
   if(discrete_quantile_residual(dist, typename Dist::value_type(0), p, c) >= 0)
      return 0;  // cdf(0) >= p already: see round_to_ceil / round_to_floor.
   //
   // What happens next depends on whether we're looking for an 
   // upper or lower quantile:
   //
   if(pp < 0.5f)
      return round_to_floor(dist, do_inverse_discrete_quantile(
         dist, 
         p, 
         c,
         (guess < 1 ? value_type(1) : (value_type)floor(guess)), 
         multiplier, 
         adder, 
         tools::equal_floor(),
         max_iter), p, c);
   // else:
   return round_to_ceil(dist, do_inverse_discrete_quantile(
      dist, 
      p, 
      c,
      (value_type)ceil(guess), 
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
      const policies::discrete_quantile<policies::integer_round_inwards>&,
      boost::math::uintmax_t& max_iter)
{
   typedef typename Dist::value_type value_type;
   BOOST_MATH_STD_USING
   typename Dist::value_type pp = c ? 1 - p : p;
   if(discrete_quantile_residual(dist, typename Dist::value_type(0), p, c) >= 0)
      return 0;  // cdf(0) >= p already: see round_to_ceil / round_to_floor.
   //
   // What happens next depends on whether we're looking for an 
   // upper or lower quantile:
   //
   if(pp < 0.5f)
      return round_to_ceil(dist, do_inverse_discrete_quantile(
         dist, 
         p, 
         c,
         ceil(guess), 
         multiplier, 
         adder, 
         tools::equal_ceil(),
         max_iter), p, c);
   // else:
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
      const policies::discrete_quantile<policies::integer_round_down>&,
      boost::math::uintmax_t& max_iter)
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
      boost::math::uintmax_t& max_iter)
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
      boost::math::uintmax_t& max_iter)
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

}}} // namespaces

#endif // BOOST_MATH_DISTRIBUTIONS_DETAIL_INV_DISCRETE_QUANTILE

