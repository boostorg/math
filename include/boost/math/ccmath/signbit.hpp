//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_SIGNBIT_HPP
#define BOOST_MATH_CCMATH_SIGNBIT_HPP

#include <cmath>
#include <limits>
#include <type_traits>
#include <boost/math/tools/is_constant_evaluated.hpp>
#include <boost/math/special_functions/detail/fp_traits.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/abs.hpp>

#if __cpp_lib_bit_cast >= 201806L
#include <bit>
#  define BOOST_MATH_BIT_CAST(T, x) std::bit_cast<T>(x)
#elif defined(__has_builtin)
#  if __has_builtin(__builtin_bit_cast)
#    define BOOST_MATH_BIT_CAST(T, x) __builtin_bit_cast(T, x)
#  endif
#endif

/*
The following error is given using Apple Clang version 13.1.6
TODO: Remove the following undef when Clang supports

ccmath_signbit_test.cpp:32:19: error: static_assert expression is not an integral constant expression
    static_assert(boost::math::ccmath::signbit(T(-1)) == true);
                  ^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
../../../boost/math/ccmath/signbit.hpp:62:24: note: constexpr bit_cast involving bit-field is not yet supported
        const auto u = BOOST_MATH_BIT_CAST(float_bits, arg);
                       ^
../../../boost/math/ccmath/signbit.hpp:20:37: note: expanded from macro 'BOOST_MATH_BIT_CAST'
#  define BOOST_MATH_BIT_CAST(T, x) __builtin_bit_cast(T, x)
                                    ^
*/

#if defined(__clang__) && defined(BOOST_MATH_BIT_CAST)
#  undef BOOST_MATH_BIT_CAST
#endif

namespace boost::math::ccmath {

namespace detail {

#ifdef BOOST_MATH_BIT_CAST

struct float_bits
{
#if BOOST_MATH_ENDIAN_LITTLE_BYTE
    unsigned mantissa : 23;
    unsigned exponent : 8;
    unsigned sign : 1;
#else // Big endian
    unsigned sign : 1;
    unsigned exponent : 8;
    unsigned mantissa : 23;
#endif 
};

template <typename T>
constexpr bool signbit_impl(T arg)
{
    if constexpr (std::is_same_v<T, float>)
    {   
        const auto u = BOOST_MATH_BIT_CAST(float_bits, arg);
        return u.sign;
    }

    if (boost::math::ccmath::isnan(arg))
    {
        return false;
    }
    
    return arg < static_cast<T>(0);
}

#else

// Typical implementations of signbit involve type punning via union and manipulating
// overflow (see libc++ or musl). Neither of these are allowed in constexpr contexts
// (technically type punning via union in general is UB in c++ but well defined in C) 
// therefore NANs and 0s are treated as positive if bit cast is unavailable

template <typename T>
constexpr bool signbit_impl(T arg)
{
    if (boost::math::ccmath::isnan(arg))
    {
        return false;
    }
    
    return arg < static_cast<T>(0);
}

#endif

}

// Return value: true if arg is negative, false if arg is 0, NAN, or positive
template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
constexpr bool signbit(Real arg)
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::detail::signbit_impl(arg);
    }
    else
    {
        using std::signbit;
        return signbit(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
constexpr bool signbit(Z arg)
{
    return boost::math::ccmath::signbit(static_cast<double>(arg));
}

} // Namespaces

#endif // BOOST_MATH_CCMATH_SIGNBIT_HPP
