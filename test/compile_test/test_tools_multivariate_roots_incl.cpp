//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header
// #includes all the files that it needs to.
//
#include <boost/math/tools/multivariate_roots.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
   using V = std::vector<double>;
   auto f = [](V const& x) { return std::make_pair(V{x[0] - 1}, V{1.0}); };
   std::uintmax_t max_iter = 10;
   check_result<V>(boost::math::tools::multivariate_newton(f, V{0.0}, 53, max_iter));
   check_result<V>(boost::math::tools::multivariate_newton(f, V{0.0}, 53));
}
