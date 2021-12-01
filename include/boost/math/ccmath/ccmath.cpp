//  (C) Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  A module containing all of the constexpr implementation of <cmath>

// Global module fragment required for non-module preprocessing
module;

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/is_constant_evaluated.hpp>

// TODO: Eventually #include for the entire STL can be replaced with import
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cinttypes>
#include <type_traits>
#include <limits>
#include <stdexcept>

export module boost.math.ccmath;

// Forward Declarations
export namespace boost::math::ccmath {

template <typename T, std::enable_if_t<!std::is_unsigned_v<T>, bool> = true>
inline constexpr T abs(T x) noexcept;

template <typename T, std::enable_if_t<std::is_unsigned_v<T>, bool> = true>
inline constexpr T abs(T x) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real ceil(Real arg) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double ceil(Z arg) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real copysign(Real mag, Real sgn) noexcept;

template <typename T1, typename T2>
inline constexpr auto copysign(T1 mag, T2 sgn) noexcept;

template <typename Z>
inline constexpr auto div(Z x, Z y) noexcept;

template <typename Z>
struct div_t
{
    Z quot;
    Z rem;
};

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real floor(Real arg) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double floor(Z arg) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real fmod(Real x, Real y) noexcept;

template <typename T1, typename T2>
inline constexpr auto fmod(T1 x, T2 y) noexcept;

template <typename T, std::enable_if_t<!std::is_integral_v<T>, bool> = true>
inline constexpr int fpclassify(T x);

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr int fpclassify(Z x);

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real frexp(Real arg, int* exp);

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double frexp(Z arg, int* exp);

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real hypot(Real x, Real y) noexcept;

template <typename T1, typename T2>
inline constexpr auto hypot(T1 x, T2 y) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr int ilogb(Real arg) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr int ilogb(Z arg) noexcept;

template <typename T>
inline constexpr bool isfinite(T x);

template <typename T>
inline constexpr bool isinf(T x);

template <typename T>
inline constexpr bool isnan(T x);

template <typename T>
inline constexpr bool isnormal(T x);

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real ldexp(Real arg, int exp) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double ldexp(Z arg, int exp) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real logb(Real arg) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double logb(Z arg) noexcept;

template <typename Real>
inline constexpr Real modf(Real x, Real* iptr);

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real remainder(Real x, Real y) noexcept;

template <typename T1, typename T2>
inline constexpr auto remainder(T1 x, T2 y) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real round(Real arg) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double round(Z arg) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real scalbln(Real arg, long exp) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double scalbln(Z arg, long exp) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real scalbn(Real arg, int exp) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double scalbn(Z arg, int exp) noexcept;

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real sqrt(Real x);

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double sqrt(Z x);

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real trunc(Real arg) noexcept;

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double trunc(Z arg) noexcept;

}

namespace boost::math::ccmath::detail {

// Detail only helper functions first
template <typename T>
inline constexpr void swap(T& x, T& y) noexcept
{
    T temp = x;
    x = y;
    y = temp;
}

template <typename T> 
inline constexpr T abs_impl(T x) noexcept
{
    return boost::math::ccmath::isnan(x) ? std::numeric_limits<T>::quiet_NaN() : 
           boost::math::ccmath::isinf(x) ? std::numeric_limits<T>::infinity() : 
           x == -0 ? T(0) :
           x == (std::numeric_limits<T>::min)() ? std::numeric_limits<T>::quiet_NaN() : 
           x > 0 ? x : -x;
}

template <typename T>
inline constexpr T ceil_impl(T arg) noexcept
{
    T result = boost::math::ccmath::floor(arg);

    if(result == arg)
    {
        return result;
    }
    else
    {
        return result + 1;
    }
}

template <typename T>
inline constexpr T copysign_impl(const T mag, const T sgn) noexcept
{
    if(sgn >= 0)
    {
        return boost::math::ccmath::abs(mag);
    }
    else
    {
        return -boost::math::ccmath::abs(mag);
    }
}

template <typename ReturnType, typename Z>
inline constexpr ReturnType div_impl(const Z x, const Z y) noexcept
{
    // std::div_t/ldiv_t/lldiv_t/imaxdiv_t can be defined as either { Z quot; Z rem; }; or { Z rem; Z quot; };
    // so don't use braced initialziation to guarantee compatibility
    ReturnType ans {0, 0};

    ans.quot = x / y;
    ans.rem = x % y;

    return ans;
}

template <typename T>
inline constexpr T floor_pos_impl(T arg) noexcept
{
    T result = 1;

    if(result < arg)
    {
        while(result < arg)
        {
            result *= 2;
        }
        while(result > arg)
        {
            --result;
        }

        return result;
    }
    else
    {
        return T(0);
    }
}

template <typename T>
inline constexpr T floor_neg_impl(T arg) noexcept
{
    T result = -1;

    if(result > arg)
    {
        while(result > arg)
        {
            result *= 2;
        }
        while(result < arg)
        {
            ++result;
        }
        if(result != arg)
        {
            --result;
        }
    }

    return result;
}

template <typename T>
inline constexpr T floor_impl(T arg) noexcept
{
    if(arg > 0)
    {
        return floor_pos_impl(arg);
    }
    else
    {
        return floor_neg_impl(arg);
    }
}

template <typename ReturnType, typename T1, typename T2>
inline constexpr ReturnType fmod_impl(T1 x, T2 y) noexcept
{
    if(x == y)
    {
        return ReturnType(0);
    }
    else
    {
        while(x >= y)
        {
            x -= y;
        }

        return static_cast<ReturnType>(x);
    }
}

template <typename Real>
inline constexpr Real frexp_zero_impl(Real arg, int* exp)
{
    *exp = 0;
    return arg;
}

template <typename Real>
inline constexpr Real frexp_impl(Real arg, int* exp)
{
    const bool negative_arg = (arg < Real(0));
    
    Real f = negative_arg ? -arg : arg;
    int e2 = 0;
    constexpr Real two_pow_32 = Real(4294967296);

    while (f >= two_pow_32)
    {
        f = f / two_pow_32;
        e2 += 32;
    }

    while(f >= Real(1))
    {
        f = f / Real(2);
        ++e2;
    }
    
    if(exp != nullptr)
    {
        *exp = e2;
    }

    return !negative_arg ? f : -f;
}

template <typename T>
inline constexpr T hypot_impl(T x, T y) noexcept
{
    x = boost::math::ccmath::abs(x);
    y = boost::math::ccmath::abs(y);

    if (y > x)
    {
        boost::math::ccmath::detail::swap(x, y);
    }

    if(x * std::numeric_limits<T>::epsilon() >= y)
    {
        return x;
    }

    T rat = y / x;
    return x * boost::math::ccmath::sqrt(1 + rat * rat);
}

template <typename Real>
inline constexpr Real ldexp_impl(Real arg, int exp) noexcept
{
    while(exp > 0)
    {
        arg *= 2;
        --exp;
    }
    while(exp < 0)
    {
        arg /= 2;
        ++exp;
    }

    return arg;
}

// The value of the exponent returned by std::logb is always 1 less than the exponent returned by 
// std::frexp because of the different normalization requirements: for the exponent e returned by std::logb, 
// |arg*r^-e| is between 1 and r (typically between 1 and 2), but for the exponent e returned by std::frexp, 
// |arg*2^-e| is between 0.5 and 1. 
template <typename T>
inline constexpr T logb_impl(T arg) noexcept
{
    int exp = 0;
    boost::math::ccmath::frexp(arg, &exp);

    return exp - 1;
}

template <typename Real>
inline constexpr Real modf_error_impl(Real x, Real* iptr)
{
    *iptr = x;
    return boost::math::ccmath::abs(x) == Real(0) ? x :
           x > Real(0) ? Real(0) : -Real(0);
}

template <typename Real>
inline constexpr Real modf_nan_impl(Real x, Real* iptr)
{
    *iptr = x;
    return x;
}

template <typename Real>
inline constexpr Real modf_impl(Real x, Real* iptr)
{
    *iptr = boost::math::ccmath::trunc(x);
    return (x - *iptr);
}

template <typename T>
inline constexpr T remainder_impl(const T x, const T y) noexcept
{
    T n = 0;
    const T fractional_part = boost::math::ccmath::modf((x / y), &n);

    if(fractional_part > T(1.0/2))
    {
        ++n;
    }
    else if(fractional_part < T(-1.0/2))
    {
        --n;
    }

    return x - n*y;
}

// Computes the nearest integer value to arg (in floating-point format), 
// rounding halfway cases away from zero, regardless of the current rounding mode.
template <typename T>
inline constexpr T round_impl(T arg) noexcept
{
    T iptr = 0;
    const T x = boost::math::ccmath::modf(arg, &iptr);
    constexpr T half = T(1)/2;

    if(x >= half && iptr > 0)
    {
        return iptr + 1;
    }
    else if(boost::math::ccmath::abs(x) >= half && iptr < 0)
    {
        return iptr - 1;
    }
    else
    {
        return iptr;
    }
}

template <typename ReturnType, typename T>
inline constexpr ReturnType int_round_impl(T arg)
{
    const T rounded_arg = round_impl(arg);

    if(rounded_arg > static_cast<T>((std::numeric_limits<ReturnType>::max)()))
    {
        if constexpr (std::is_same_v<ReturnType, long long>)
        {
            throw std::domain_error("Rounded value cannot be represented by a long long type without overflow");
        }
        else
        {
            throw std::domain_error("Rounded value cannot be represented by a long type without overflow");
        }
    }
    else
    {
        return static_cast<ReturnType>(rounded_arg);
    }
}

template <typename Real, typename Z>
inline constexpr Real scalbn_impl(Real arg, Z exp) noexcept
{
    while(exp > 0)
    {
        arg *= FLT_RADIX;
        --exp;
    }
    while(exp < 0)
    {
        arg /= FLT_RADIX;
        ++exp;
    }

    return arg;
}

template <typename Real>
inline constexpr Real sqrt_impl_2(Real x, Real s, Real s2)
{
    return !(s < s2) ? s2 : sqrt_impl_2(x, (x / s + s) / 2, s);
}

template <typename Real>
inline constexpr Real sqrt_impl_1(Real x, Real s)
{
    return sqrt_impl_2(x, (x / s + s) / 2, s);
}

template <typename Real>
inline constexpr Real sqrt_impl(Real x)
{
    return sqrt_impl_1(x, x > 1 ? x : Real(1));
}

template <typename T>
inline constexpr T trunc_impl(T arg) noexcept
{
    return (arg > 0) ? boost::math::ccmath::floor(arg) : boost::math::ccmath::ceil(arg);
}

} // Namespace boost::math::ccmath::detail

// Useable Functions

export namespace boost::math::ccmath {

template <typename T, std::enable_if_t<!std::is_unsigned_v<T>, bool> = true>
inline constexpr T abs(T x) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return detail::abs_impl<T>(x);
    }
    else
    {
        using std::abs;
        return abs(x);
    }
}

// If abs() is called with an argument of type X for which is_unsigned_v<X> is true and if X
// cannot be converted to int by integral promotion (7.3.7), the program is ill-formed.
template <typename T, std::enable_if_t<std::is_unsigned_v<T>, bool> = true>
inline constexpr T abs(T x) noexcept
{
    if constexpr (std::is_convertible_v<T, int>)
    {
        return detail::abs_impl<int>(static_cast<int>(x));
    }
    else
    {
        static_assert(sizeof(T) == 0, "Taking the absolute value of an unsigned value not covertible to int is UB.");
        return T(0); // Unreachable, but suppresses warnings
    }
}

inline constexpr long int labs(long int j) noexcept
{
    return boost::math::ccmath::abs(j);
}

inline constexpr long long int llabs(long long int j) noexcept
{
    return boost::math::ccmath::abs(j);
}

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real ceil(Real arg) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               boost::math::ccmath::detail::ceil_impl(arg);
    }
    else
    {
        using std::ceil;
        return ceil(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double ceil(Z arg) noexcept
{
    return boost::math::ccmath::ceil(static_cast<double>(arg));
}

inline constexpr float ceilf(float arg) noexcept
{
    return boost::math::ccmath::ceil(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double ceill(long double arg) noexcept
{
    return boost::math::ccmath::ceil(arg);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real copysign(Real mag, Real sgn) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(mag))
    {
        return boost::math::ccmath::detail::copysign_impl(mag, sgn);
    }
    else
    {
        using std::copysign;
        return copysign(mag, sgn);
    }
}

template <typename T1, typename T2>
inline constexpr auto copysign(T1 mag, T2 sgn) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(mag))
    {
        // If the type is an integer (e.g. epsilon == 0) then set the epsilon value to 1 so that type is at a minimum 
        // cast to double
        constexpr auto T1p = std::numeric_limits<T1>::epsilon() > 0 ? std::numeric_limits<T1>::epsilon() : 1;
        constexpr auto T2p = std::numeric_limits<T2>::epsilon() > 0 ? std::numeric_limits<T2>::epsilon() : 1;
        
        using promoted_type = 
                              #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                              std::conditional_t<T1p <= LDBL_EPSILON && T1p <= T2p, T1,
                              std::conditional_t<T2p <= LDBL_EPSILON && T2p <= T1p, T2,
                              #endif
                              std::conditional_t<T1p <= DBL_EPSILON && T1p <= T2p, T1,
                              std::conditional_t<T2p <= DBL_EPSILON && T2p <= T1p, T2, double
                              #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                              >>>>;
                              #else
                              >>;
                              #endif

        return boost::math::ccmath::copysign(promoted_type(mag), promoted_type(sgn));
    }
    else
    {
        using std::copysign;
        return copysign(mag, sgn);
    }
}

inline constexpr float copysignf(float mag, float sgn) noexcept
{
    return boost::math::ccmath::copysign(mag, sgn);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double copysignl(long double mag, long double sgn) noexcept
{
    return boost::math::ccmath::copysign(mag, sgn);
}
#endif

template <typename Z>
inline constexpr auto div(Z x, Z y) noexcept
{
    if constexpr (std::is_same_v<Z, int>)
    {
        return detail::div_impl<std::div_t>(x, y);
    }
    else if constexpr (std::is_same_v<Z, long>)
    {
        return detail::div_impl<std::ldiv_t>(x, y);
    }
    else if constexpr (std::is_same_v<Z, long long>)
    {
        return detail::div_impl<std::lldiv_t>(x, y);
    }
    else if constexpr (std::is_same_v<Z, std::intmax_t>)
    {
        return detail::div_impl<std::imaxdiv_t>(x, y);
    }
    else
    {
        return detail::div_impl<boost::math::ccmath::div_t<Z>>(x, y);
    }
}

inline constexpr std::ldiv_t ldiv(long x, long y) noexcept
{
    return detail::div_impl<std::ldiv_t>(x, y);
}

inline constexpr std::lldiv_t lldiv(long long x, long long y) noexcept
{
    return detail::div_impl<std::lldiv_t>(x, y);
}

inline constexpr std::imaxdiv_t imaxdiv(std::intmax_t x, std::intmax_t y) noexcept
{
    return detail::div_impl<std::imaxdiv_t>(x, y);
}

template <typename T>
inline constexpr auto fabs(T x) noexcept
{
    return boost::math::ccmath::abs(x);
}

inline constexpr float fabsf(float x) noexcept
{
    return boost::math::ccmath::abs(x);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double fabsl(long double x) noexcept
{
    return boost::math::ccmath::abs(x);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real floor(Real arg) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               boost::math::ccmath::detail::floor_impl(arg);
    }
    else
    {
        using std::floor;
        return floor(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double floor(Z arg) noexcept
{
    return boost::math::ccmath::floor(static_cast<double>(arg));
}

inline constexpr float floorf(float arg) noexcept
{
    return boost::math::ccmath::floor(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double floorl(long double arg) noexcept
{
    return boost::math::ccmath::floor(arg);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real fmod(Real x, Real y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return boost::math::ccmath::abs(x) == Real(0) && y != Real(0) ? x :
               boost::math::ccmath::isinf(x) && !boost::math::ccmath::isnan(y) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::abs(y) == Real(0) && !boost::math::ccmath::isnan(x) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::isinf(y) && boost::math::ccmath::isfinite(x) ? x :
               boost::math::ccmath::isnan(x) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::isnan(y) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::detail::fmod_impl<Real>(x, y);
    }
    else
    {
        using std::fmod;
        return fmod(x, y);
    }
}

template <typename T1, typename T2>
inline constexpr auto fmod(T1 x, T2 y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        // If the type is an integer (e.g. epsilon == 0) then set the epsilon value to 1 so that type is at a minimum 
        // cast to double
        constexpr auto T1p = std::numeric_limits<T1>::epsilon() > 0 ? std::numeric_limits<T1>::epsilon() : 1;
        constexpr auto T2p = std::numeric_limits<T2>::epsilon() > 0 ? std::numeric_limits<T2>::epsilon() : 1;
        
        using promoted_type = 
                              #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                              std::conditional_t<T1p <= LDBL_EPSILON && T1p <= T2p, T1,
                              std::conditional_t<T2p <= LDBL_EPSILON && T2p <= T1p, T2,
                              #endif
                              std::conditional_t<T1p <= DBL_EPSILON && T1p <= T2p, T1,
                              std::conditional_t<T2p <= DBL_EPSILON && T2p <= T1p, T2, double
                              #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                              >>>>;
                              #else
                              >>;
                              #endif

        return boost::math::ccmath::fmod(promoted_type(x), promoted_type(y));
    }
    else
    {
        using std::fmod;
        return fmod(x, y);
    }
}

inline constexpr float fmodf(float x, float y) noexcept
{
    return boost::math::ccmath::fmod(x, y);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double fmodl(long double x, long double y) noexcept
{
    return boost::math::ccmath::fmod(x, y);
}
#endif

template <typename T, std::enable_if_t<!std::is_integral_v<T>, bool> = true>
inline constexpr int fpclassify(T x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return boost::math::ccmath::isnan(x) ? FP_NAN :
               boost::math::ccmath::isinf(x) ? FP_INFINITE :
               boost::math::ccmath::abs(x) == T(0) ? FP_ZERO :
               boost::math::ccmath::abs(x) > 0 && boost::math::ccmath::abs(x) < (std::numeric_limits<T>::min)() ? FP_SUBNORMAL : FP_NORMAL;
    }
    else
    {
        using std::fpclassify;
        return fpclassify(x);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr int fpclassify(Z x)
{
    return boost::math::ccmath::fpclassify(static_cast<double>(x));
}

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real frexp(Real arg, int* exp)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return arg == Real(0)  ? detail::frexp_zero_impl(arg, exp) : 
               arg == Real(-0) ? detail::frexp_zero_impl(arg, exp) :
               boost::math::ccmath::isinf(arg) ? detail::frexp_zero_impl(arg, exp) : 
               boost::math::ccmath::isnan(arg) ? detail::frexp_zero_impl(arg, exp) :
               boost::math::ccmath::detail::frexp_impl(arg, exp);
    }
    else
    {
        using std::frexp;
        return frexp(arg, exp);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double frexp(Z arg, int* exp)
{
    return boost::math::ccmath::frexp(static_cast<double>(arg), exp);
}

inline constexpr float frexpf(float arg, int* exp)
{
    return boost::math::ccmath::frexp(arg, exp);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double frexpl(long double arg, int* exp)
{
    return boost::math::ccmath::frexp(arg, exp);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real hypot(Real x, Real y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return boost::math::ccmath::abs(x) == Real(0) ? boost::math::ccmath::abs(y) :
               boost::math::ccmath::abs(y) == Real(0) ? boost::math::ccmath::abs(x) :
               boost::math::ccmath::isinf(x) ? std::numeric_limits<Real>::infinity() :
               boost::math::ccmath::isinf(y) ? std::numeric_limits<Real>::infinity() :
               boost::math::ccmath::isnan(x) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::isnan(y) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::detail::hypot_impl(x, y);
    }
    else
    {
        using std::hypot;
        return hypot(x, y);
    }
}

template <typename T1, typename T2>
inline constexpr auto hypot(T1 x, T2 y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        // If the type is an integer (e.g. epsilon == 0) then set the epsilon value to 1 so that type is at a minimum 
        // cast to double
        constexpr auto T1p = std::numeric_limits<T1>::epsilon() > 0 ? std::numeric_limits<T1>::epsilon() : 1;
        constexpr auto T2p = std::numeric_limits<T2>::epsilon() > 0 ? std::numeric_limits<T2>::epsilon() : 1;
        
        using promoted_type = 
                              #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                              std::conditional_t<T1p <= LDBL_EPSILON && T1p <= T2p, T1,
                              std::conditional_t<T2p <= LDBL_EPSILON && T2p <= T1p, T2,
                              #endif
                              std::conditional_t<T1p <= DBL_EPSILON && T1p <= T2p, T1,
                              std::conditional_t<T2p <= DBL_EPSILON && T2p <= T1p, T2, double
                              #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                              >>>>;
                              #else
                              >>;
                              #endif

        return boost::math::ccmath::hypot(promoted_type(x), promoted_type(y));
    }
    else
    {
        using std::hypot;
        return hypot(x, y);
    }
}

inline constexpr float hypotf(float x, float y) noexcept
{
    return boost::math::ccmath::hypot(x, y);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double hypotl(long double x, long double y) noexcept
{
    return boost::math::ccmath::hypot(x, y);
}
#endif

// If arg is not zero, infinite, or NaN, the value returned is exactly equivalent to static_cast<int>(std::logb(arg))
template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr int ilogb(Real arg) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? FP_ILOGB0 :
               boost::math::ccmath::isinf(arg) ? INT_MAX :
               boost::math::ccmath::isnan(arg) ? FP_ILOGBNAN :
               static_cast<int>(boost::math::ccmath::logb(arg));
    }
    else
    {
        using std::ilogb;
        return ilogb(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr int ilogb(Z arg) noexcept
{
    return boost::math::ccmath::ilogb(static_cast<double>(arg));
}

inline constexpr int ilogbf(float arg) noexcept
{
    return boost::math::ccmath::ilogb(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr int ilogbl(long double arg) noexcept
{
    return boost::math::ccmath::ilogb(arg);
}
#endif

template <typename T>
inline constexpr bool isfinite(T x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        // bool isfinite (IntegralType arg) is a set of overloads accepting the arg argument of any integral type
        // equivalent to casting the integral argument arg to double (e.g. static_cast<double>(arg))
        if constexpr (std::is_integral_v<T>)
        {
            return !boost::math::ccmath::isinf(static_cast<double>(x)) && !boost::math::ccmath::isnan(static_cast<double>(x));
        }
        else
        {
            return !boost::math::ccmath::isinf(x) && !boost::math::ccmath::isnan(x);
        }
    }
    else
    {
        using std::isfinite;

        if constexpr (!std::is_integral_v<T>)
        {
            return isfinite(x);
        }
        else
        {
            return isfinite(static_cast<double>(x));
        }
    }
}

template <typename T>
inline constexpr bool isinf(T x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return x == std::numeric_limits<T>::infinity() || -x == std::numeric_limits<T>::infinity();
    }
    else
    {
        using std::isinf;
        
        if constexpr (!std::is_integral_v<T>)
        {
            return isinf(x);
        }
        else
        {
            return isinf(static_cast<double>(x));
        }
    }
}

template <typename T>
inline constexpr bool isnan(T x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return x != x;
    }
    else
    {
        using std::isnan;

        if constexpr (!std::is_integral_v<T>)
        {
            return isnan(x);
        }
        else
        {
            return isnan(static_cast<double>(x));
        }
    }
}

template <typename T>
inline constexpr bool isnormal(T x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {   
        return x == T(0) ? false :
               boost::math::ccmath::isinf(x) ? false :
               boost::math::ccmath::isnan(x) ? false :
               boost::math::ccmath::abs(x) < (std::numeric_limits<T>::min)() ? false : true;
    }
    else
    {
        using std::isnormal;

        if constexpr (!std::is_integral_v<T>)
        {
            return isnormal(x);
        }
        else
        {
            return isnormal(static_cast<double>(x));
        }
    }
}

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real ldexp(Real arg, int exp) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               boost::math::ccmath::detail::ldexp_impl(arg, exp);
    }
    else
    {
        using std::ldexp;
        return ldexp(arg, exp);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double ldexp(Z arg, int exp) noexcept
{
    return boost::math::ccmath::ldexp(static_cast<double>(arg), exp);
}

inline constexpr float ldexpf(float arg, int exp) noexcept
{
    return boost::math::ccmath::ldexp(arg, exp);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double ldexpl(long double arg, int exp) noexcept
{
    return boost::math::ccmath::ldexp(arg, exp);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real logb(Real arg) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? -std::numeric_limits<Real>::infinity() :
               boost::math::ccmath::isinf(arg) ? std::numeric_limits<Real>::infinity() :
               boost::math::ccmath::isnan(arg) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::detail::logb_impl(arg);
    }
    else
    {
        using std::logb;
        return logb(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double logb(Z arg) noexcept
{
    return boost::math::ccmath::logb(static_cast<double>(arg));
}

inline constexpr float logbf(float arg) noexcept
{
    return boost::math::ccmath::logb(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double logbl(long double arg) noexcept
{
    return boost::math::ccmath::logb(arg);
}
#endif

template <typename Real>
inline constexpr Real modf(Real x, Real* iptr)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return boost::math::ccmath::abs(x) == Real(0) ? detail::modf_error_impl(x, iptr) :
               boost::math::ccmath::isinf(x) ? detail::modf_error_impl(x, iptr) :
               boost::math::ccmath::isnan(x) ? detail::modf_nan_impl(x, iptr) :
               boost::math::ccmath::detail::modf_impl(x, iptr);
    }
    else
    {
        using std::modf;
        return modf(x, iptr);
    }
}

inline constexpr float modff(float x, float* iptr)
{
    return boost::math::ccmath::modf(x, iptr);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double modfl(long double x, long double* iptr)
{
    return boost::math::ccmath::modf(x, iptr);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real remainder(Real x, Real y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return boost::math::ccmath::isinf(x) && !boost::math::ccmath::isnan(y) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::abs(y) == Real(0) && !boost::math::ccmath::isnan(x) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::isnan(x) || boost::math::ccmath::isnan(y) ? std::numeric_limits<Real>::quiet_NaN() :
               boost::math::ccmath::detail::remainder_impl<Real>(x, y);
    }
    else
    {
        using std::remainder;
        return remainder(x, y);
    }
}

template <typename T1, typename T2>
inline constexpr auto remainder(T1 x, T2 y) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        // If the type is an integer (e.g. epsilon == 0) then set the epsilon value to 1 so that type is at a minimum 
        // cast to double
        constexpr auto T1p = std::numeric_limits<T1>::epsilon() > 0 ? std::numeric_limits<T1>::epsilon() : 1;
        constexpr auto T2p = std::numeric_limits<T2>::epsilon() > 0 ? std::numeric_limits<T2>::epsilon() : 1;
        
        using promoted_type = 
                              #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                              std::conditional_t<T1p <= LDBL_EPSILON && T1p <= T2p, T1,
                              std::conditional_t<T2p <= LDBL_EPSILON && T2p <= T1p, T2,
                              #endif
                              std::conditional_t<T1p <= DBL_EPSILON && T1p <= T2p, T1,
                              std::conditional_t<T2p <= DBL_EPSILON && T2p <= T1p, T2, double
                              #ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
                              >>>>;
                              #else
                              >>;
                              #endif

        return boost::math::ccmath::remainder(promoted_type(x), promoted_type(y));
    }
    else
    {
        using std::remainder;
        return remainder(x, y);
    }
}

inline constexpr float remainderf(float x, float y) noexcept
{
    return boost::math::ccmath::remainder(x, y);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double remainderl(long double x, long double y) noexcept
{
    return boost::math::ccmath::remainder(x, y);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real round(Real arg) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               boost::math::ccmath::detail::round_impl(arg);
    }
    else
    {
        using std::round;
        return round(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double round(Z arg) noexcept
{
    return boost::math::ccmath::round(static_cast<double>(arg));
}

inline constexpr float roundf(float arg) noexcept
{
    return boost::math::ccmath::round(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double roundl(long double arg) noexcept
{
    return boost::math::ccmath::round(arg);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr long lround(Real arg)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? 0l :
               boost::math::ccmath::isinf(arg) ? 0l :
               boost::math::ccmath::isnan(arg) ? 0l :
               boost::math::ccmath::detail::int_round_impl<long>(arg);
    }
    else
    {
        using std::lround;
        return lround(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr long lround(Z arg)
{
    return boost::math::ccmath::lround(static_cast<double>(arg));
}

inline constexpr long lroundf(float arg)
{
    return boost::math::ccmath::lround(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long lroundl(long double arg)
{
    return boost::math::ccmath::lround(arg);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr long long llround(Real arg)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? 0ll :
               boost::math::ccmath::isinf(arg) ? 0ll :
               boost::math::ccmath::isnan(arg) ? 0ll :
               boost::math::ccmath::detail::int_round_impl<long long>(arg);
    }
    else
    {
        using std::llround;
        return llround(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr long llround(Z arg)
{
    return boost::math::ccmath::llround(static_cast<double>(arg));
}

inline constexpr long long llroundf(float arg)
{
    return boost::math::ccmath::llround(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long long llroundl(long double arg)
{
    return boost::math::ccmath::llround(arg);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real scalbln(Real arg, long exp) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               boost::math::ccmath::detail::scalbn_impl(arg, exp);
    }
    else
    {
        using std::scalbln;
        return scalbln(arg, exp);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double scalbln(Z arg, long exp) noexcept
{
    return boost::math::ccmath::scalbln(static_cast<double>(arg), exp);
}

inline constexpr float scalblnf(float arg, long exp) noexcept
{
    return boost::math::ccmath::scalbln(arg, exp);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double scalblnl(long double arg, long exp) noexcept
{
    return boost::math::ccmath::scalbln(arg, exp);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real scalbn(Real arg, int exp) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               boost::math::ccmath::detail::scalbn_impl(arg, exp);
    }
    else
    {
        using std::scalbn;
        return scalbn(arg, exp);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double scalbn(Z arg, int exp) noexcept
{
    return boost::math::ccmath::scalbn(static_cast<double>(arg), exp);
}

inline constexpr float scalbnf(float arg, int exp) noexcept
{
    return boost::math::ccmath::scalbn(arg, exp);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double scalbnl(long double arg, int exp) noexcept
{
    return boost::math::ccmath::scalbn(arg, exp);
}
#endif

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real sqrt(Real x)
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(x))
    {
        return boost::math::ccmath::isnan(x) ? std::numeric_limits<Real>::quiet_NaN() : 
               boost::math::ccmath::isinf(x) ? std::numeric_limits<Real>::infinity() : 
               detail::sqrt_impl<Real>(x);
    }
    else
    {
        using std::sqrt;
        return sqrt(x);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double sqrt(Z x)
{
    return detail::sqrt_impl<double>(static_cast<double>(x));
}

template <typename Real, std::enable_if_t<!std::is_integral_v<Real>, bool> = true>
inline constexpr Real trunc(Real arg) noexcept
{
    if(BOOST_MATH_IS_CONSTANT_EVALUATED(arg))
    {
        return boost::math::ccmath::abs(arg) == Real(0) ? arg :
               boost::math::ccmath::isinf(arg) ? arg :
               boost::math::ccmath::isnan(arg) ? arg :
               boost::math::ccmath::detail::trunc_impl(arg);
    }
    else
    {
        using std::trunc;
        return trunc(arg);
    }
}

template <typename Z, std::enable_if_t<std::is_integral_v<Z>, bool> = true>
inline constexpr double trunc(Z arg) noexcept
{
    return boost::math::ccmath::trunc(static_cast<double>(arg));
}

inline constexpr float truncf(float arg) noexcept
{
    return boost::math::ccmath::trunc(arg);
}

#ifndef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
inline constexpr long double truncl(long double arg) noexcept
{
    return boost::math::ccmath::trunc(arg);
}
#endif

}
