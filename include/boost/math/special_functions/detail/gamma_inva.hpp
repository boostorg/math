//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// This is not a complete header file, it is included by gamma.hpp
// after it has defined it's definitions.  This inverts the incomplete
// gamma functions P and Q on the first parameter "a" using a generic
// root finding algorithm (TOMS Algorithm 748).
//

#ifndef BOOST_MATH_SP_DETAIL_GAMMA_INVA
#define BOOST_MATH_SP_DETAIL_GAMMA_INVA

#include <boost/math/tools/toms748_solve.hpp>
#include <boost/cstdint.hpp>

namespace boost{ namespace math{ namespace detail{

template <class T>
struct gamma_inva_t
{
   gamma_inva_t(T z_, T p_, bool invert_) : z(z_), p(p_), invert(invert_) {}
   T operator()(T a)
   {
      return invert ? p - boost::math::gamma_q(a, z) : boost::math::gamma_p(a, z) - p;
   }
private:
   T z, p;
   bool invert;
};

template <class T>
T gamma_inva_imp(const T& z, const T& p, const T& q)
{
   using namespace std;  // for ADL of std lib math functions
   //
   // Special cases first:
   //
   if(p == 0)
   {
      return tools::max_value<T>();
   }
   if(q == 0)
   {
      return tools::min_value<T>();
   }
   //
   // Function object, this is the functor whose root
   // we have to solve:
   //
   gamma_inva_t<T> f(z, (p < q) ? p : q, (p < q) ? false : true);
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
   if(z > 1.1)
   {
      guess = z;
   }
   else if(z > 0.5)
   {
      guess = z * 1.2f;
   }
   else
   {
      guess = -0.4f / log(z);
   }
   //
   // Max iterations permitted:
   //
   boost::uintmax_t max_iter = 200;
   std::pair<T, T> r = bracket_and_solve_root(f, guess, static_cast<T>(2), false, tol, max_iter);
   if(max_iter >= 200)
      tools::logic_error<T>(BOOST_CURRENT_FUNCTION, "Unable to locate the root within a reasonable number of iterations, closest approximation so far was %1%", r.first);
   return (r.first + r.second) / 2;
}

} // namespace detail

template <class T>
inline T gamma_p_inva(T x, T p)
{
   return detail::gamma_inva_imp(x, p, 1 - p);
}

template <class T>
inline T gamma_q_inva(T x, T q)
{
   return detail::gamma_inva_imp(x, 1 - q, q);
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SP_DETAIL_GAMMA_INVA

