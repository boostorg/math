//  Copyright John Maddock 2005-2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_PRECISION_INCLUDED
#define BOOST_MATH_TOOLS_PRECISION_INCLUDED

#include <boost/limits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>

#include <iostream>
#include <iomanip>

namespace boost{ namespace math
{
namespace tools
{

template <class T>
inline int digits(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   return std::numeric_limits<T>::digits;
}

template <class T>
inline T max_value(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   return (std::numeric_limits<T>::max)();
}

template <class T>
inline T min_value(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
#endif
   return (std::numeric_limits<T>::min)();
}

namespace detail{
//
// Logarithmic limits come next, note that although
// we can compute these from the log of the max value
// that is not in general thread safe (if we cache the value)
// so it's better to specialise these:
//
// For type float first:
//
template <class T>
inline T log_max_value(const mpl::int_<128>& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
{
   return 88.0f;
}

template <class T>
inline T log_min_value(const mpl::int_<128>& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
{
   return -87.0f;
}
//
// Now double:
//
template <class T>
inline T log_max_value(const mpl::int_<1024>& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
{
   return 709.0;
}

template <class T>
inline T log_min_value(const mpl::int_<1024>& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
{
   return -708.0;
}
//
// 80 and 128-bit long doubles:
//
template <class T>
inline T log_max_value(const mpl::int_<16384>& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
{
   return 11356.0L;
}

template <class T>
inline T log_min_value(const mpl::int_<16384>& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
{
   return -11355.0L;
}

template <class T, class U>
inline T log_max_value(const U& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
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

template <class T, class U>
inline T log_min_value(const U& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
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
inline T epsilon(const mpl::true_& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
{
   return std::numeric_limits<T>::epsilon();
}

template <class T>
inline T epsilon(const mpl::false_& BOOST_APPEND_EXPLICIT_TEMPLATE_TYPE(T))
{
   using namespace std;  // for ADL of std names
   static const T eps = ldexp(static_cast<T>(1), 1-tools::digits<T>());
   return eps;
}

} // namespace detail

template <class T>
inline T log_max_value(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
   return detail::log_max_value<T>(mpl::int_<std::numeric_limits<T>::max_exponent>());
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
   using namespace std;
   static const T val = log((std::numeric_limits<T>::max)());
   return val;
#endif
}

template <class T>
inline T log_min_value(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   BOOST_STATIC_ASSERT( ::std::numeric_limits<T>::is_specialized);
   return detail::log_min_value<T>(mpl::int_<std::numeric_limits<T>::max_exponent>());
#else
   BOOST_ASSERT(::std::numeric_limits<T>::is_specialized);
   using namespace std;
   static const T val = log((std::numeric_limits<T>::min)());
   return val;
#endif
}

template <class T>
inline T epsilon(BOOST_EXPLICIT_TEMPLATE_TYPE(T))
{
#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
   return detail::epsilon<T>(mpl::bool_< ::std::numeric_limits<T>::is_specialized>());
#else
   return ::std::numeric_limits<T>::is_specialized ?
      detail::epsilon<T>(mpl::true_()) :
      detail::epsilon<T>(mpl::false_());
#endif
}

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_PRECISION_INCLUDED
