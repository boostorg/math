//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/beta.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/beta.hpp>

template float boost::math::beta<float>(float, float);
template double boost::math::beta<double>(double, double);
template long double boost::math::beta<long double>(long double, long double);

template float boost::math::ibeta<float>(float, float, float);
template double boost::math::ibeta<double>(double, double, double);
template long double boost::math::ibeta<long double>(long double, long double, long double);

template float boost::math::ibeta_inv<float>(float, float, float);
template double boost::math::ibeta_inv<double>(double, double, double);
template long double boost::math::ibeta_inv<long double>(long double, long double, long double);

template float boost::math::ibeta_inva<float>(float, float, float);
template double boost::math::ibeta_inva<double>(double, double, double);
template long double boost::math::ibeta_inva<long double>(long double, long double, long double);

