//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/erf.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/erf.hpp>

template float boost::math::erf<float>(float);
template double boost::math::erf<double>(double);
template long double boost::math::erf<long double>(long double);

template float boost::math::erfc<float>(float);
template double boost::math::erfc<double>(double);
template long double boost::math::erfc<long double>(long double);

template float boost::math::erf_inv<float>(float);
template double boost::math::erf_inv<double>(double);
template long double boost::math::erf_inv<long double>(long double);

template float boost::math::erfc_inv<float>(float);
template double boost::math::erfc_inv<double>(double);
template long double boost::math::erfc_inv<long double>(long double);

