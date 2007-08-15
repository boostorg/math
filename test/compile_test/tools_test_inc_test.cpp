//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/tools/test.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/tools/test.hpp>
#include <boost/array.hpp>

template float boost::math::tools::relative_error<float>(float a, float b);

#define A boost::array<boost::array<double, 2>, 2>
typedef double (*F1)(const boost::array<double, 2>&);
typedef F1 F2;

template boost::math::tools::test_result<
   boost::math::tools::calculate_result_type<A>::value_type> 
      boost::math::tools::test<A, F1, F2>(const A& a, F1 test_func, F2 expect_func);

