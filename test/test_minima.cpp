//  Copyright John Maddock 2006.
//  Copyright Paul A. Bristow 2007.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/tools/minima.hpp>
#include <boost/test/included/test_exec_monitor.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/gamma.hpp>

template <class T>
struct poly_test
{
   // minima is at (3,4):
   T operator()(const T& v)
   {
      T a = v - 3;
      return 3 * a * a + 4;
   }
};

template <class T>
void test_minima(T, const char* /* name */)
{
   std::pair<T, T> m = boost::math::tools::brent_find_minima(poly_test<T>(), T(-10), T(10), 50);
   BOOST_CHECK_CLOSE(m.first, T(3), T(0.001));
   BOOST_CHECK_CLOSE(m.second, T(4), T(0.001));

   T (*fp)(T);
   fp = boost::math::lgamma;

   m = boost::math::tools::brent_find_minima(fp, T(0.5), T(10), 50);
   BOOST_CHECK_CLOSE(m.first, T(1.461632), T(0.1));
   fp = boost::math::tgamma;
   m = boost::math::tools::brent_find_minima(fp, T(0.5), T(10), 50);
   BOOST_CHECK_CLOSE(m.first, T(1.461632), T(0.1));
}

int test_main(int, char* [])
{
   test_minima(0.1f, "float");
   test_minima(0.1, "double");
   test_minima(0.1L, "long double");
   return 0;
}


