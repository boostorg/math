
//  Copyright (c) 2011 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_BIG_CONSTANT_HPP
#define BOOST_MATH_TOOLS_BIG_CONSTANT_HPP

#include <boost/math/tools/config.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/type_traits/is_convertible.hpp>

namespace boost{ namespace math{ namespace tools{

template <class T>
inline T make_big_value(long double v, const char*, mpl::true_ const&, mpl::false_ const&)
{
   return static_cast<T>(v);
}
template <class T>
inline T make_big_value(long double v, const char*, mpl::true_ const&, mpl::true_ const&)
{
   return static_cast<T>(v);
}
template <class T>
inline T make_big_value(long double, const char* s, mpl::false_ const&, mpl::false_ const&)
{
   return boost::lexical_cast<T>(s);
}
template <class T>
inline const char* make_big_value(long double, const char* s, mpl::false_ const&, mpl::true_ const&)
{
   return s;
}

#define BOOST_MATH_BIG_CONSTANT(T, D, x)\
   boost::math::tools::make_big_value<T>(BOOST_JOIN(x, L), BOOST_STRINGIZE(x), mpl::bool_<D <= std::numeric_limits<long double>::digits>(), boost::is_convertible<const char*, T>())

}}} // namespaces

#endif

