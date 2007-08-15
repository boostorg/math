//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/fpclassify.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/fpclassify.hpp>

template int boost::math::fpclassify<float>(float);
template int boost::math::fpclassify<double>(double);
template int boost::math::fpclassify<long double>(long double);

template bool boost::math::isfinite<float>(float);
template bool boost::math::isfinite<double>(double);
template bool boost::math::isfinite<long double>(long double);

template bool boost::math::isinf<float>(float);
template bool boost::math::isinf<double>(double);
template bool boost::math::isinf<long double>(long double);

template bool boost::math::isnormal<float>(float);
template bool boost::math::isnormal<double>(double);
template bool boost::math::isnormal<long double>(long double);



