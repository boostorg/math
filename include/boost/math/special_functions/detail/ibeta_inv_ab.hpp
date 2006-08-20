//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// This is not a complete header file, it is included by beta.hpp
// after it has defined it's definitions.  This inverts the incomplete
// beta functions ibeta and ibetac on the first parameters "a"
// and "b" using a generic root finding algorithm (TOMS Algorithm 748).
//

#ifndef BOOST_MATH_SP_DETAIL_BETA_INV_AB
#define BOOST_MATH_SP_DETAIL_BETA_INV_AB

#include <boost/math/tools/toms748_solve.hpp>
#include <boost/cstdint.hpp>

namespace boost{ namespace math{ namespace detail{

template <class T>
struct beta_inv_ab_t
{
   beta_inv_ab_t(T b_, T z_, T p_, bool invert_, bool swap_ab_) : b(b_), z(z_), p(p_), invert(invert_), swap_ab(swap_ab_) {}
   T operator()(T a)
   {
      return invert ? 
         p - boost::math::ibetac(swap_ab ? b : a, swap_ab ? a : b, z) 
         : boost::math::ibeta(swap_ab ? b : a, swap_ab ? a : b, z) - p;
   }
private:
   T b, z, p;
   bool invert, swap_ab;
};

template <class T>
T ibeta_inv_ab_imp(const T& b, const T& z, const T& p, const T& q, bool swap_ab)
{
   using namespace std;  // for ADL of std lib math functions
   //
   // Special cases first:
   //
   if(p == 0)
   {
      return swap_ab ? tools::min_value<T>() : tools::max_value<T>();
   }
   if(q == 0)
   {
      return swap_ab ? tools::max_value<T>() : tools::min_value<T>();
   }
   //
   // Function object, this is the functor whose root
   // we have to solve:
   //
   beta_inv_ab_t<T> f(b, z, (p < q) ? p : q, (p < q) ? false : true, swap_ab);
   //
   // Tolerance: full precision.
   //
   tools::eps_tolerance<T> tol(tools::digits<T>());
   //
   // Now figure out a starting guess for what a may be, 
   // we'll start out with a value that'll put p or q
   // right bang in the middle of their range, the functions
   // are quite sensitive so we should need too many steps
   // to bracket the root from there:
   //
   T guess;
   if((p < q) != swap_ab)
   {
      guess = b * 2;
   }
   else
   {
      guess = b / 2;
   }
   //
   // Max iterations permitted:
   //
   boost::uintmax_t max_iter = 200;
   std::pair<T, T> r = bracket_and_solve_root(f, guess, static_cast<T>(5), swap_ab ? true : false, tol, max_iter);
   if(max_iter >= 200)
      tools::logic_error<T>(BOOST_CURRENT_FUNCTION, "Unable to locate the root within a reasonable number of iterations, closest approximation so far was %1%", r.first);
   return (r.first + r.second) / 2;
}

} // namespace detail

template <class T>
inline T ibeta_inva(T b, T x, T p)
{
   return detail::ibeta_inv_ab_imp(b, x, p, 1 - p, false);
}

template <class T>
inline T ibetac_inva(T b, T x, T q)
{
   return detail::ibeta_inv_ab_imp(b, x, 1 - q, q, false);
}

template <class T>
inline T ibeta_invb(T b, T x, T p)
{
   return detail::ibeta_inv_ab_imp(b, x, p, 1 - p, true);
}

template <class T>
inline T ibetac_invb(T b, T x, T q)
{
   return detail::ibeta_inv_ab_imp(b, x, 1 - q, q, true);
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SP_DETAIL_BETA_INV_AB

