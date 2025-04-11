//  (C) Copyright Christopher Kormanyos 2025.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/core/lightweight_test.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <cmath>
#include <limits>

auto main() -> int
{
  using float_type = float;

  const float_type x { 0x1.03ebbp-128F };

  const float_type ctrl { 0x1.41085ep+127F };

  const float_type result = ::boost::math::cyl_neumann(-1, x);

  const float_type tol = std::numeric_limits<float_type>::epsilon() * 16;

  using std::fabs;

  BOOST_TEST(result > 0);
  BOOST_TEST(fabs(1 - (result / ctrl)) < tol);

  return boost::report_errors();
}
