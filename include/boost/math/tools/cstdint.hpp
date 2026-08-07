//  Copyright (c) 2024 Matt Borland
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_CSTDINT
#define BOOST_MATH_TOOLS_CSTDINT

#include <boost/math/tools/config.hpp>


#ifdef BOOST_MATH_ENABLE_CUDA

#include <cuda/std/cstdint>

namespace boost {
namespace math {

BOOST_MATH_EXPORT using cuda::std::int8_t;
BOOST_MATH_EXPORT using cuda::std::int16_t;
BOOST_MATH_EXPORT using cuda::std::int32_t;
BOOST_MATH_EXPORT using cuda::std::int64_t;

BOOST_MATH_EXPORT using cuda::std::int_fast8_t;
BOOST_MATH_EXPORT using cuda::std::int_fast16_t;
BOOST_MATH_EXPORT using cuda::std::int_fast32_t;
BOOST_MATH_EXPORT using cuda::std::int_fast64_t;

BOOST_MATH_EXPORT using cuda::std::int_least8_t;
BOOST_MATH_EXPORT using cuda::std::int_least16_t;
BOOST_MATH_EXPORT using cuda::std::int_least32_t;
BOOST_MATH_EXPORT using cuda::std::int_least64_t;

BOOST_MATH_EXPORT using cuda::std::intmax_t;
BOOST_MATH_EXPORT using cuda::std::intptr_t;

BOOST_MATH_EXPORT using cuda::std::uint8_t;
BOOST_MATH_EXPORT using cuda::std::uint16_t;
BOOST_MATH_EXPORT using cuda::std::uint32_t;
BOOST_MATH_EXPORT using cuda::std::uint64_t;

BOOST_MATH_EXPORT using cuda::std::uint_fast8_t;
BOOST_MATH_EXPORT using cuda::std::uint_fast16_t;
BOOST_MATH_EXPORT using cuda::std::uint_fast32_t;
BOOST_MATH_EXPORT using cuda::std::uint_fast64_t;

BOOST_MATH_EXPORT using cuda::std::uint_least8_t;
BOOST_MATH_EXPORT using cuda::std::uint_least16_t;
BOOST_MATH_EXPORT using cuda::std::uint_least32_t;
BOOST_MATH_EXPORT using cuda::std::uint_least64_t;

BOOST_MATH_EXPORT using cuda::std::uintmax_t;
BOOST_MATH_EXPORT using cuda::std::uintptr_t;

using size_t = unsigned long;

#else

#ifndef BOOST_MATH_BUILD_MODULE
#include <cstdint>
#endif

namespace boost {
namespace math {

BOOST_MATH_EXPORT using std::int8_t;
BOOST_MATH_EXPORT using std::int16_t;
BOOST_MATH_EXPORT using std::int32_t;
BOOST_MATH_EXPORT using std::int64_t;

BOOST_MATH_EXPORT using std::int_fast8_t;
BOOST_MATH_EXPORT using std::int_fast16_t;
BOOST_MATH_EXPORT using std::int_fast32_t;
BOOST_MATH_EXPORT using std::int_fast64_t;

BOOST_MATH_EXPORT using std::int_least8_t;
BOOST_MATH_EXPORT using std::int_least16_t;
BOOST_MATH_EXPORT using std::int_least32_t;
BOOST_MATH_EXPORT using std::int_least64_t;

BOOST_MATH_EXPORT using std::intmax_t;
BOOST_MATH_EXPORT using std::intptr_t;

BOOST_MATH_EXPORT using std::uint8_t;
BOOST_MATH_EXPORT using std::uint16_t;
BOOST_MATH_EXPORT using std::uint32_t;
BOOST_MATH_EXPORT using std::uint64_t;

BOOST_MATH_EXPORT using std::uint_fast8_t;
BOOST_MATH_EXPORT using std::uint_fast16_t;
BOOST_MATH_EXPORT using std::uint_fast32_t;
BOOST_MATH_EXPORT using std::uint_fast64_t;

BOOST_MATH_EXPORT using std::uint_least8_t;
BOOST_MATH_EXPORT using std::uint_least16_t;
BOOST_MATH_EXPORT using std::uint_least32_t;
BOOST_MATH_EXPORT using std::uint_least64_t;

BOOST_MATH_EXPORT using std::uintmax_t;
BOOST_MATH_EXPORT using std::uintptr_t;

BOOST_MATH_EXPORT using std::size_t;

#endif

} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_CSTDINT
