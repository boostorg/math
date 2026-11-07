//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/fresnel.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/fresnel.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void compile_and_link_test()
{
   check_result<float>(boost::math::fresnel_c(f));
   check_result<double>(boost::math::fresnel_c(d));
   check_result<double>(boost::math::fresnel_s(d));
   check_result<double>(boost::math::fresnel_s(d, boost::math::policies::policy<>()));
   check_result<double>(boost::math::fresnel_cis(d).real());
   check_result<float>(boost::math::fresnel_cis(f, boost::math::policies::policy<>()).imag());
}
