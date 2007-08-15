//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/erf.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/bessel.hpp>

template float boost::math::cyl_bessel_j<float>(float, float);
template double boost::math::cyl_bessel_j<double>(double, double);
template long double boost::math::cyl_bessel_j<long double>(long double, long double);

template float boost::math::cyl_neumann<float>(float, float);
template double boost::math::cyl_neumann<double>(double, double);
template long double boost::math::cyl_neumann<long double>(long double, long double);

template float boost::math::cyl_bessel_i<float>(float, float);
template double boost::math::cyl_bessel_i<double>(double, double);
template long double boost::math::cyl_bessel_i<long double>(long double, long double);

template float boost::math::cyl_bessel_k<float>(float, float);
template double boost::math::cyl_bessel_k<double>(double, double);
template long double boost::math::cyl_bessel_k<long double>(long double, long double);

template float boost::math::sph_bessel<float>(unsigned, float);
template double boost::math::sph_bessel<double>(unsigned, double);
template long double boost::math::sph_bessel<long double>(unsigned, long double);

template float boost::math::sph_neumann<float>(unsigned, float);
template double boost::math::sph_neumann<double>(unsigned, double);
template long double boost::math::sph_neumann<long double>(unsigned, long double);

