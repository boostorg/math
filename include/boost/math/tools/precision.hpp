//  Copyright John Maddock 2005-2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_PRECISION_INCLUDED
#define BOOST_MATH_TOOLS_PRECISION_INCLUDED

#include <boost/limits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <iostream>
#include <iomanip>

namespace boost{ namespace math
{
namespace tools
{

template <class T>
int digits(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   return std::numeric_limits<T>::digits;
}

template <class T>
T max_value(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   return (std::numeric_limits<T>::max)();
}

template <class T>
T min_value(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   return (std::numeric_limits<T>::min)();
}

template <class T>
T log_max_value(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   using namespace std;
   static const T val = log((std::numeric_limits<T>::max)());
   return val;
}

template <class T>
T log_min_value(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   using namespace std;
   static const T val = log((std::numeric_limits<T>::min)());
   return val;
}

template <class T>
T epsilon(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   return std::numeric_limits<T>::epsilon();
}

template <class T>
void setprecision(std::ostream& os, T, int p)
{
   os << std::setprecision(p);
}

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_PRECISION_INCLUDED
