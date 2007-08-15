//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/spherical_harmonic.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/spherical_harmonic.hpp>

template std::complex<float> boost::math::spherical_harmonic<float>(unsigned, int, float, float);
template std::complex<double> boost::math::spherical_harmonic<double>(unsigned, int, double, double);
template std::complex<long double> boost::math::spherical_harmonic<long double>(unsigned, int, long double, long double);

template float boost::math::spherical_harmonic_r<float>(unsigned, int, float, float);
template double boost::math::spherical_harmonic_r<double>(unsigned, int, double, double);
template long double boost::math::spherical_harmonic_r<long double>(unsigned, int, long double, long double);

template float boost::math::spherical_harmonic_i<float>(unsigned, int, float, float);
template double boost::math::spherical_harmonic_i<double>(unsigned, int, double, double);
template long double boost::math::spherical_harmonic_i<long double>(unsigned, int, long double, long double);



