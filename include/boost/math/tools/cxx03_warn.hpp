//  Copyright (c) 2020 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_CXX03_WARN_HPP
#define BOOST_MATH_TOOLS_CXX03_WARN_HPP

#include <boost/config/pragma_message.hpp>

#if defined(BOOST_NO_CXX11_NOEXCEPT) || defined(BOOST_NO_CXX11_RVALUE_REFERENCES) || defined(BOOST_NO_SFINAE_EXPR)\
      || defined(BOOST_NO_CXX11_AUTO_DECLARATIONS) || defined(BOOST_NO_CXX11_LAMBDAS) || defined(BOOST_NO_CXX11_UNIFIED_INITIALIZATION_SYNTAX)\
      || defined(BOOST_NO_CXX11_HDR_TUPLE) || defined(BOOST_NO_CXX11_HDR_INITIALIZER_LIST) || defined(BOOST_NO_CXX11_HDR_CHRONO)\
      || defined(BOOST_NO_CXX11_THREAD_LOCAL) || defined(BOOST_NO_CXX11_CONSTEXPR) || defined(BOOST_NO_CXX11_NULLPTR)\
      || defined(BOOST_NO_CXX11_NUMERIC_LIMITS) || defined(BOOST_NO_CXX11_DECLTYPE) || defined(BOOST_NO_CXX11_HDR_ARRAY)\
      || defined(BOOST_NO_CXX11_ALLOCATOR) || defined(BOOST_NO_CXX11_HDR_ATOMIC) || defined(BOOST_NO_CXX11_NOEXCEPT)\
      || defined(BOOST_NO_CXX11_EXPLICIT_CONVERSION_OPERATORS)
//
// The above list includes everything we use, plus a few we're likely to use soon.
// As from March 2020, C++03 support is deprecated, and as from March 2021 will be removed,
// so mark up as such:
//
#if (defined(_MSC_VER) || defined(__GNUC__)) && !defined(BOOST_MATH_DISABLE_DEPRECATED_03_WARNING)
BOOST_PRAGMA_MESSAGE("CAUTION: One or more C++11 features were found to be unavailable")
BOOST_PRAGMA_MESSAGE("CAUTION: Compiling Boost.Math in pre-C++11 conformance modes is now deprecated and will be removed from March 2021.")
BOOST_PRAGMA_MESSAGE("CAUTION: Define BOOST_MATH_DISABLE_DEPRECATED_03_WARNING to suppress this message.")
#endif
#endif

#endif
