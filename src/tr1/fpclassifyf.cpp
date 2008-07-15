//  Copyright John Maddock 2008.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#define BOOST_MATH_TR1_SOURCE
#include <boost/math/tr1.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/math/special_functions/sign.hpp>
#include "c_policy.hpp"

namespace boost{ namespace math{ namespace tr1{

#define NO_MACRO_EXPAND /**/

template<> bool BOOST_MATH_TR1_DECL signbit<float> NO_MACRO_EXPAND(float x)
{
   return static_cast<bool>((boost::math::signbit)(x)); 
}

template<> int BOOST_MATH_TR1_DECL fpclassify<float> NO_MACRO_EXPAND(float x)
{
   return (boost::math::fpclassify)(x);
}

template<> bool BOOST_MATH_TR1_DECL isfinite<float> NO_MACRO_EXPAND(float x)
{
   return (boost::math::isfinite)(x);
}

template<> bool BOOST_MATH_TR1_DECL isinf<float> NO_MACRO_EXPAND(float x)
{
   return (boost::math::isinf)(x);
}

template<> bool BOOST_MATH_TR1_DECL isnan<float> NO_MACRO_EXPAND(float x)
{
   return (boost::math::isnan)(x);
}

template<> bool BOOST_MATH_TR1_DECL isnormal<float> NO_MACRO_EXPAND(float x)
{
   return (boost::math::isnormal)(x);
}

}}}
