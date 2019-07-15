//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_CARDINAL_B_SPLINE_HPP
#define BOOST_MATH_SPECIAL_CARDINAL_B_SPLINE_HPP

#include <array>
#include <cmath>

namespace boost { namespace math {

namespace detail {

  template<class Real>
  inline Real B1(Real x)
  {
    if (x < 0)
    {
      return B1(-x);
    }
    if (x < Real(1))
    {
      return 1 - x;
    }
    return Real(0);
  }
}

template<unsigned n, typename Real>
Real cardinal_b_spline(Real x) {
  if (x < 0) {
      // All B-splines are even functions:
      return cardinal_b_spline<n>(-x);
  }

  // Someday we'll be able to 'if constexpr' on many compilers.
  // Not today though! In any case, it's not a big performance hit.
  if /*constexpr*/ (n==0)
  {
    if (x < Real(1)/Real(2)) {
      return Real(1);
    }
    else if (x == Real(1)/Real(2)) {
      return Real(1)/Real(2);
    }
    else {
      return Real(0);
    }
  }

  if /*constexpr*/ (n==1) {
    return detail::B1(x);
  }

  Real supp_max = (n+1)/Real(2);
  if (x >= supp_max)
  {
    return Real(0);
  }

  // Fill v with values of B1:
  // At most two of these terms are nonzero, and at least 1.
  // There is only one non-zero term when n is odd and x = 0.
  std::array<Real, n> v;
  Real z = x + 1 - supp_max;
  for (unsigned i = 0; i < n; ++i)
  {
      v[i] = detail::B1(z);
      z += 1;
  }

  Real smx = supp_max - x;
  for (unsigned j = 2; j <= n; ++j)
  {
    Real a = (j + 1 - smx);
    Real b = smx;
    for(unsigned k = 0; k <= n - j; ++k)
    {
      v[k] = (a*v[k+1] + b*v[k])/Real(j);
      a += 1;
      b -= 1;
    }
  }

  return v[0];
}

template<unsigned n, class Real>
Real forward_cardinal_b_spline(Real x) {
    return cardinal_b_spline<n>(x - (n+1)/Real(2));
}

}}
#endif
