//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Regular use of std::array functions can not be used on 
//  GPU platforms like CUDA since they are missing the __device__ marker
//  Alias as needed to get correct support

#ifndef BOOST_MATH_TOOLS_ARRAY_HPP
#define BOOST_MATH_TOOLS_ARRAY_HPP

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_ENABLE_CUDA

#include <cuda/std/array>

BOOST_MATH_NAMESPACE_BEGIN

using cuda::std::array;

BOOST_MATH_NAMESPACE_END

#else

#ifndef BOOST_MATH_BUILD_MODULE
#include <array>
#endif

BOOST_MATH_NAMESPACE_BEGIN

using std::array;

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_ENABLE_CUDA

#endif // BOOST_MATH_TOOLS_ARRAY_HPP
