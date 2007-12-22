//  Copyright John Maddock 2007.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_MODF_HPP
#define BOOST_MATH_MODF_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions/trunc.hpp>

namespace boost{ namespace math{

template <class T>
inline T modf(const T& v, T* ipart)
{
   *ipart = trunc(v);
   return v - *ipart;
}

template <class T>
inline T modf(const T& v, int* ipart)
{
   *ipart = itrunc(v);
   return v - *ipart;
}

template <class T>
inline T modf(const T& v, long* ipart)
{
   *ipart = ltrunc(v);
   return v - *ipart;
}

#ifdef BOOST_HAS_LONG_LONG
template <class T>
inline T modf(const T& v, long long* ipart)
{
   *ipart = lltrunc(v);
   return v - *ipart;
}
#endif

}} // namespaces

#endif // BOOST_MATH_MODF_HPP
