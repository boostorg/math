//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/special_functions/laguerre.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/special_functions/laguerre.hpp>

template float boost::math::laguerre<float>(unsigned, float);
template double boost::math::laguerre<double>(unsigned, double);
template long double boost::math::laguerre<long double>(unsigned, long double);

typedef boost::math::policies::policy<> def_pol;

template float boost::math::laguerre<float, def_pol>(unsigned, unsigned, float, const def_pol&);
template double boost::math::laguerre<double, def_pol>(unsigned, unsigned, double, const def_pol&);
template long double boost::math::laguerre<long double, def_pol>(unsigned, unsigned, long double, const def_pol&);

