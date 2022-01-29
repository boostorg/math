//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Macros that substitute for STL concepts or typename depending on availability of <concepts>

#ifndef BOOST_MATH_TOOLS_CONCEPTS_HPP
#define BOOST_MATH_TOOLS_CONCEPTS_HPP

#if (__cplusplus > 202000L || _MSVC_LANG > 202000L)
#  if __has_include(<concepts>)
#    include <concepts>
#    define BOOST_MATH_FLOATING_POINT_TYPE std::floating_point
#  else
#    define BOOST_MATH_FLOATING_POINT_TYPE typename
#  endif
#else
#  define BOOST_MATH_FLOATING_POINT_TYPE typename
#endif

#endif // BOOST_MATH_TOOLS_CONCEPTS_HPP
