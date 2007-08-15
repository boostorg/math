//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/distributions/normal.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/distributions/normal.hpp>

template class boost::math::normal_distribution<float, boost::math::policies::policy<> >;
template class boost::math::normal_distribution<double, boost::math::policies::policy<> >;
template class boost::math::normal_distribution<long double, boost::math::policies::policy<> >;
