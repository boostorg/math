//  Copyright John Maddock 2007.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TRUNC_HPP
#define BOOST_MATH_TRUNC_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/math/tools/config.hpp>

namespace boost{ namespace math{

template <class T>
inline T trunc(const T& v)
{
   BOOST_MATH_STD_USING
   return (v >= 0) ? floor(v) : ceil(v);
}
//
// The following functions will not compile unless T has an
// implicit convertion to the integer types.  For user-defined
// number types this will likely not be the case.  In that case
// these functions should either be specialized for the UDT in
// question, or else overloads should be placed in the same 
// namespace as the UDT: these will then be found via argument
// dependent lookup.  See our concept archetypes for examples.
//
template <class T>
inline int itrunc(const T& v)
{
   return static_cast<int>(boost::math::trunc(v));
}

template <class T>
inline long ltrunc(const T& v)
{
   return static_cast<long int>(boost::math::trunc(v));
}

#ifdef BOOST_HAS_LONG_LONG

template <class T>
inline long long lltrunc(const T& v)
{
   return static_cast<long long>(boost::math::trunc(v));
}

#endif

}} // namespaces

#endif // BOOST_MATH_TRUNC_HPP
