//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Regular use of <type_traits> is not compatible with CUDA
//  Adds aliases to unify the support
//  Also adds convience overloads like is_same_v so we don't have to wait for C++17

#ifndef BOOST_MATH_TOOLS_TYPE_TRAITS
#define BOOST_MATH_TOOLS_TYPE_TRAITS

#include <boost/math/tools/config.hpp>

namespace boost {
namespace math {

#ifdef BOOST_MATH_CUDA_ENABLED

#include <cuda/std/type_traits>

using cuda::std::is_void;
using cuda::std::is_integral;
using cuda::std::enable_if_t;

#else // STD versions

#include <type_traits>

using std::is_void;
using std::is_integral;
using std::enable_if_t;

#endif 

template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_void_v = boost::math::is_void<T>::value;

template <typename T>
BOOST_MATH_INLINE_CONSTEXPR bool is_integral_v = boost::math::is_integral<T>::value;

} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_TYPE_TRAITS
