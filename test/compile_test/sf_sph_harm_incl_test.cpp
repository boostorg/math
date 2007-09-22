//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/spherical_harmonic.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/spherical_harmonic.hpp>
//
// Note this header includes no other headers, this is
// important if this test is to be meaningful:
//
#include "test_compile_result.hpp"

void check()
{
   check_result<std::complex<float> >(boost::math::spherical_harmonic<float>(u, i, f, f));
   check_result<std::complex<double> >(boost::math::spherical_harmonic<double>(u, i, d, d));
   check_result<std::complex<long double> >(boost::math::spherical_harmonic<long double>(u, i, l, l));

   check_result<float>(boost::math::spherical_harmonic_r<float>(u, i, f, f));
   check_result<double>(boost::math::spherical_harmonic_r<double>(u, i, d, d));
   check_result<long double>(boost::math::spherical_harmonic_r<long double>(u, i, l, l));

   check_result<float>(boost::math::spherical_harmonic_i<float>(u, i, f, f));
   check_result<double>(boost::math::spherical_harmonic_i<double>(u, i, d, d));
   check_result<long double>(boost::math::spherical_harmonic_i<long double>(u, i, l, l));
}


