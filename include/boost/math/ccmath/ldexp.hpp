//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_LDEXP_HPP
#define BOOST_MATH_CCMATH_LDEXP_HPP

#include <cmath>
#include <limits>
#include <type_traits>
#include <boost/math/tools/is_constant_evaluated.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/ccmath/abs.hpp>
#include <boost/math/ccmath/isinf.hpp>
#include <boost/math/ccmath/isnan.hpp>
#include <boost/math/ccmath/isfinite.hpp>

namespace boost::math::ccmath {

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real ldexp(Real arg, int exp)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               exp == 0 ? arg :
               arg * boost::math::pow<exp>(2); // Todo replace with boost::math::ccmath::exp2
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
