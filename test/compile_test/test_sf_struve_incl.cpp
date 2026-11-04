//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/struve.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/struve.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
   check_result<float>(boost::math::struve_h(f, f));
   check_result<double>(boost::math::struve_h(d, d));
   check_result<double>(boost::math::struve_l(d, d));
   double derivative;
   check_result<double>(boost::math::struve_h(d, d, &derivative));
   check_result<float>(boost::math::struve_l(f, f, static_cast<float*>(nullptr), boost::math::policies::policy<>()));
   check_result<double>(boost::math::struve_h(d, d, boost::math::policies::policy<>()));
   check_result<double>(boost::math::struve_l(d, d, boost::math::policies::policy<>()));
}
