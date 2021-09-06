//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_LDEXP_HPP
#define BOOST_MATH_CCMATH_LDEXP_HPP

#include <cmath>
#include <limits>
#include <type_traits>
#include <stdexcept>
#include <boost/math/tools/is_constant_evaluated.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/ccmath/abs.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/isfinite.hpp>

namespace boost::math::ccmath {

namespace detail {

template <typename T>
inline constexpr T negative_exp(int exp)
{
    int abs_exp = boost::math::ccmath::abs(exp);

    // std::feclearexcept and std::fetestexcept are not constexpr
    // Simple test if approaching underflow
    constexpr T denorm_min_limit = std::numeric_limits<T>::denorm_min() * 2;
    T result = 2;

    while(abs_exp >= 0)
    {
        if(result < denorm_min_limit)
        {
            throw std::underflow_error("Exponent leads to a value that is smaller than minimum positive subnormal value");
        }
        
        result /= 2;
        --abs_exp;
    }

    return result;
}

template <typename T>
inline constexpr T positive_exp(int exp)
{
    constexpr T max_limit = (std::numeric_limits<T>::max)() / 2;
    T result = 2;

    while((exp-1) > 0)
    {
        if(result > max_limit)
        {
            throw std::overflow_error("Exponent leads to a value that is larger than the maximum allowable value without overflow");
        }

        result *= 2;
        --exp;
    }

    return result;
}

// A rudimentary implementation of exp2 that only deals with int exponents for simplicity
template <typename T>
inline constexpr T intexp2(int exp)
{
    if(exp > 0)
    {
        return positive_exp<T>(exp);
    }
    else
    {
        return negative_exp<T>(exp);
    }
}

} // Namespace detail

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real ldexp(Real arg, int exp)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               exp == 0 ? arg :
               arg * boost::math::ccmath::detail::intexp2<Real>(exp);
    }
    else
    {
        using std::ldexp;
        return ldexp(arg, exp);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double ldexp(Z arg, int exp)
{
    return boost::math::ccmath::ldexp(static_cast<double>(arg), exp);
}

inline constexpr float ldexpf(float arg, int exp)
{
    return boost::math::ccmath::ldexp(arg, exp);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double ldexpl(long double arg, int exp)
{
    return boost::math::ccmath::ldexp(arg, exp);
}
#endif

} // Namespaces

#endif // BOOST_MATH_CCMATH_LDEXP_HPP
