//  (C) Copyright John Maddock 2010.
//  (C) Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TUPLE_HPP_INCLUDED
#define BOOST_MATH_TUPLE_HPP_INCLUDED

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_ENABLE_CUDA

#include <thrust/pair.h>
#include <thrust/tuple.h>

namespace boost { 
namespace math {

using thrust::pair;
using thrust::tuple;

using thrust::make_pair;
using thrust::make_tuple;

using thrust::tie;
using thrust::get;

using thrust::tuple_size;
using thrust::tuple_element;

} // namespace math
} // namespace boost

#else

#include <tuple>

namespace boost { 
namespace math {

using ::std::tuple;
using ::std::pair;

// [6.1.3.2] Tuple creation functions
using ::std::ignore;
using ::std::make_tuple;
using ::std::tie;
using ::std::get;

// [6.1.3.3] Tuple helper classes
using ::std::tuple_size;
using ::std::tuple_element;

// Pair helpers
using ::std::make_pair;

} // namespace math
} // namespace boost

#endif // BOOST_MATH_ENABLE_CUDA

#endif
