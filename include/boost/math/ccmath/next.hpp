//  (C) Copyright John Maddock 2008 - 2022.
//  (C) Copyright Matt Borland 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_CCMATH_NEXT_HPP
#define BOOST_MATH_CCMATH_NEXT_HPP

#include <cmath>
#include <cfloat>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <stdexcept>
#include <boost/math/tools/assert.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/is_constant_evaluated.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/traits.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/math/ccmath/ilogb.hpp>
#include <boost/math/ccmath/ldexp.hpp>
#include <boost/math/ccmath/scalbln.hpp>
#include <boost/math/ccmath/round.hpp>
#include <boost/math/ccmath/fabs.hpp>
#include <boost/math/ccmath/fpclassify.hpp>
#include <boost/math/ccmath/isfinite.hpp>
#include <boost/math/ccmath/fmod.hpp>
#include <boost/math/ccmath/detail/minmax.hpp>

namespace boost::math::concepts {
    class real_concept;
    class std_real_concept;
}

namespace boost::math::ccmath {

namespace detail {

template <typename T>
struct has_hidden_guard_digits;
template <>
struct has_hidden_guard_digits<float> : public std::false_type {};
template <>
struct has_hidden_guard_digits<double> : public std::false_type {};
template <>
struct has_hidden_guard_digits<long double> : public std::false_type {};
#ifdef BOOST_HAS_FLOAT128
template <>
struct has_hidden_guard_digits<__float128> : public std::false_type {};
#endif
template <>
struct has_hidden_guard_digits<boost::math::concepts::real_concept> : public std::false_type {};
template <>
struct has_hidden_guard_digits<boost::math::concepts::std_real_concept> : public std::false_type {};

template <typename T, bool b>
struct has_hidden_guard_digits_10 : public std::false_type {};
template <typename T>
struct has_hidden_guard_digits_10<T, true> : public std::integral_constant<bool, (std::numeric_limits<T>::digits10 != std::numeric_limits<T>::max_digits10)> {};

template <typename T>
struct has_hidden_guard_digits 
    : public has_hidden_guard_digits_10<T, 
    std::numeric_limits<T>::is_specialized
    && (std::numeric_limits<T>::radix == 10) >
{};

template <typename T>
constexpr T& normalize_value(const T& val, const std::false_type&) { return val; }
template <typename T>
constexpr T normalize_value(const T& val, const std::true_type&) 
{
    static_assert(std::numeric_limits<T>::is_specialized, "Type T must be specialized.");
    static_assert(std::numeric_limits<T>::radix != 2, "Type T must be specialized.");

    std::intmax_t shift = static_cast<std::intmax_t>(std::numeric_limits<T>::digits) - static_cast<std::intmax_t>(boost::math::ccmath::ilogb(val)) - 1;
    T result = boost::math::ccmath::scalbn(val, shift);
    result = boost::math::ccmath::round(result);
    return boost::math::ccmath::scalbn(result, -shift); 
}

template <typename T>
constexpr T get_smallest_value(const std::true_type&)
{
    //
    // numeric_limits lies about denorms being present - particularly
    // when this can be turned on or off at runtime, as is the case
    // when using the SSE2 registers in DAZ or FTZ mode.
    //
    static constexpr T m = std::numeric_limits<T>::denorm_min();
    return ((tools::min_value<T>() / 2) == 0) ? tools::min_value<T>() : m;
}

template <typename T>
constexpr T get_smallest_value(const std::false_type&)
{
    return tools::min_value<T>();
}

template <typename T>
constexpr T get_smallest_value()
{
    return get_smallest_value<T>(std::integral_constant<bool, std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::has_denorm == std::denorm_present)>());
}

//
// Returns the smallest value that won't generate denorms when
// we calculate the value of the least-significant-bit:
//
template <typename T>
constexpr T get_min_shift_value();

template <typename T>
struct min_shift_initializer
{
    struct init
    {
        init()
        {
            do_init();
        }
        static constexpr void do_init()
        {
            get_min_shift_value<T>();
        }
        constexpr void force_instantiate() const {}
    };

    static constexpr init initializer;
    static constexpr void force_instantiate()
    {
        initializer.force_instantiate();
    }
};

template <typename T>
constexpr typename min_shift_initializer<T>::init min_shift_initializer<T>::initializer;

template <typename T>
constexpr T calc_min_shifted(const std::true_type&)
{
   return boost::math::ccmath::ldexp(tools::min_value<T>(), tools::digits<T>() + 1);
}

template <typename T>
constexpr T calc_min_shifted(const std::false_type&)
{
   static_assert(std::numeric_limits<T>::is_specialized, "Type T must be specialized.");
   static_assert(std::numeric_limits<T>::radix != 2, "Type T must be specialized.");

   return boost::math::ccmath::scalbn(tools::min_value<T>(), std::numeric_limits<T>::digits + 1);
}


template <typename T>
constexpr T get_min_shift_value()
{
   static const T val = calc_min_shifted<T>(std::integral_constant<bool, !std::numeric_limits<T>::is_specialized || std::numeric_limits<T>::radix == 2>());
   min_shift_initializer<T>::force_instantiate();

   return val;
}

template <typename T, bool b = boost::math::tools::detail::has_backend_type_v<T>>
struct exponent_type
{
    using type = int;
};

template <typename T>
struct exponent_type<T, true>
{
    using type = typename T::backend_type::exponent_type;
};

template <typename T, bool b = boost::math::tools::detail::has_backend_type_v<T>>
using exponent_type_t = typename exponent_type<T>::type;

template <typename T, typename Policy>
constexpr T float_next_imp(const T& val, const std::true_type&, const Policy& pol)
{
    using exponent_type = exponent_type_t<T>;
    
    exponent_type expon {};
    static const char* function = "float_next<%1%>(%1%)";

    int fpclass = boost::math::ccmath::fpclassify(val);

    if ((fpclass == FP_NAN) || (fpclass == FP_INFINITE))
    {
        if (val < 0)
        {
            return -tools::max_value<T>();
        }
        return policies::raise_domain_error<T>(
            function,
            "Argument must be finite, but got %1%", val, pol);
    }

    if (val >= tools::max_value<T>())
    {
        return policies::raise_overflow_error<T>(function, 0, pol);
    }

    if (val == 0)
    {
        return detail::get_smallest_value<T>();
    }

    if ((fpclass != FP_SUBNORMAL) && (fpclass != FP_ZERO) 
        && (boost::math::ccmath::fabs(val) < detail::get_min_shift_value<T>()) 
        && (val != -tools::min_value<T>()))
    {
        //
        // Special case: if the value of the least significant bit is a denorm, and the result
        // would not be a denorm, then shift the input, increment, and shift back.
        // This avoids issues with the Intel SSE2 registers when the FTZ or DAZ flags are set.
        //
        return boost::math::ccmath::ldexp(float_next(static_cast<T>(boost::math::ccmath::ldexp(val, 2 * tools::digits<T>())), pol), -2 * tools::digits<T>());
    }

    if (-0.5f == boost::math::ccmath::frexp(val, &expon))
    {
        --expon; // reduce exponent when val is a power of two, and negative.
    }
    T diff = boost::math::ccmath::ldexp(static_cast<T>(1), expon - tools::digits<T>());
    if(diff == 0)
    {
        diff = detail::get_smallest_value<T>();
    }
    return val + diff;
}

//
// Special version for some base other than 2:
//
template <typename T, typename Policy>
constexpr T float_next_imp(const T& val, const std::false_type&, const Policy& pol)
{
    using exponent_type = exponent_type_t<T>;

    static_assert(std::numeric_limits<T>::is_specialized, "Type T must be specialized.");
    static_assert(std::numeric_limits<T>::radix != 2, "Type T must be specialized.");

    exponent_type expon {};
    static const char* function = "float_next<%1%>(%1%)";

    int fpclass = boost::math::ccmath::fpclassify(val);

    if ((fpclass == FP_NAN) || (fpclass == FP_INFINITE))
    {
        if (val < 0)
        {
            return -tools::max_value<T>();
        }

        return policies::raise_domain_error<T>(
            function,
            "Argument must be finite, but got %1%", val, pol);
    }

    if(val >= tools::max_value<T>())
    {
        return policies::raise_overflow_error<T>(function, 0, pol);
    }

    if(val == 0)
    {
        return detail::get_smallest_value<T>();
    }

    if ((fpclass != FP_SUBNORMAL) && (fpclass != FP_ZERO) 
        && (boost::math::ccmath::fabs(val) < detail::get_min_shift_value<T>()) 
        && (val != -tools::min_value<T>()))
    {
        //
        // Special case: if the value of the least significant bit is a denorm, and the result
        // would not be a denorm, then shift the input, increment, and shift back.
        // This avoids issues with the Intel SSE2 registers when the FTZ or DAZ flags are set.
        //
        return boost::math::ccmath::scalbn(float_next(static_cast<T>(boost::math::ccmath::scalbn(val, 2 * std::numeric_limits<T>::digits)), pol), -2 * std::numeric_limits<T>::digits);
    }

    expon = 1 + boost::math::ccmath::ilogb(val);
    if(-1 == boost::math::ccmath::scalbn(val, -expon) * std::numeric_limits<T>::radix)
    {
        --expon; // reduce exponent when val is a power of base, and negative.
    }

    T diff = boost::math::ccmath::scalbn(static_cast<T>(1), expon - std::numeric_limits<T>::digits);
    if(diff == 0)
    {
        diff = detail::get_smallest_value<T>();
    }

    return val + diff;
}

} // namespace detail

template <typename T, typename Policy>
constexpr tools::promote_args_t<T> float_next(const T& val, const Policy& pol)
{
    using result_type = tools::promote_args_t<T>;
    return detail::float_next_imp(detail::normalize_value(static_cast<result_type>(val), typename detail::has_hidden_guard_digits<result_type>::type()), std::integral_constant<bool, !std::numeric_limits<result_type>::is_specialized || (std::numeric_limits<result_type>::radix == 2)>(), pol);
}

template <typename T>
constexpr tools::promote_args_t<T> float_next(const T& val)
{
    return float_next(val, policies::policy<>());
}

namespace detail {

template <typename T, typename Policy>
constexpr T float_prior_imp(const T& val, const std::true_type&, const Policy& pol)
{
    using exponent_type = exponent_type_t<T>;

    exponent_type expon {};
    static const char* function = "float_prior<%1%>(%1%)";

    int fpclass = boost::math::ccmath::fpclassify(val);

    if ((fpclass == FP_NAN) || (fpclass == FP_INFINITE))
    {
        if (val > 0)
        {
            return tools::max_value<T>();
        }

        return policies::raise_domain_error<T>(
            function,
            "Argument must be finite, but got %1%", val, pol);
    }

    if (val <= -tools::max_value<T>())
    {
        return -policies::raise_overflow_error<T>(function, 0, pol);
    }

    if (val == 0)
    {
        return -detail::get_smallest_value<T>();
    }

    if ((fpclass != FP_SUBNORMAL) && (fpclass != FP_ZERO) 
        && (boost::math::ccmath::fabs(val) < detail::get_min_shift_value<T>()) 
        && (val != tools::min_value<T>()))
    {
        //
        // Special case: if the value of the least significant bit is a denorm, and the result
        // would not be a denorm, then shift the input, increment, and shift back.
        // This avoids issues with the Intel SSE2 registers when the FTZ or DAZ flags are set.
        //
        return boost::math::ccmath::ldexp(float_prior(static_cast<T>(boost::math::ccmath::ldexp(val, 2 * tools::digits<T>())), pol), -2 * tools::digits<T>());
    }

    if(T remain = boost::math::ccmath::frexp(val, &expon); remain == 0.5f)
    {
        --expon; // when val is a power of two we must reduce the exponent
    }

    T diff = boost::math::ccmath::ldexp(static_cast<T>(1), expon - tools::digits<T>());
    if(diff == 0)
    {
        diff = detail::get_smallest_value<T>();
    }

    return val - diff;
}

//
// Special version for bases other than 2:
//
template <typename T, typename Policy>
constexpr T float_prior_imp(const T& val, const std::false_type&, const Policy& pol)
{
    using exponent_type = exponent_type_t<T>;

    static_assert(std::numeric_limits<T>::is_specialized, "Type T must be specialized.");
    static_assert(std::numeric_limits<T>::radix != 2, "Type T must be specialized.");

    exponent_type expon {};
    static const char* function = "float_prior<%1%>(%1%)";

    int fpclass = boost::math::ccmath::fpclassify(val);

    if ((fpclass == FP_NAN) || (fpclass == FP_INFINITE))
    {
        if (val > 0)
        {
            return tools::max_value<T>();
        }

        return policies::raise_domain_error<T>(
            function,
            "Argument must be finite, but got %1%", val, pol);
    }

    if (val <= -tools::max_value<T>())
    {
        return -policies::raise_overflow_error<T>(function, 0, pol);
    }

    if (val == 0)
    {
        return -detail::get_smallest_value<T>();
    }

    if ((fpclass != FP_SUBNORMAL) && (fpclass != FP_ZERO) 
        && (boost::math::ccmath::fabs(val) < detail::get_min_shift_value<T>()) 
        && (val != tools::min_value<T>()))
    {
        //
        // Special case: if the value of the least significant bit is a denorm, and the result
        // would not be a denorm, then shift the input, increment, and shift back.
        // This avoids issues with the Intel SSE2 registers when the FTZ or DAZ flags are set.
        //
        return boost::math::ccmath::scalbn(float_prior(static_cast<T>(boost::math::ccmath::scalbn(val, 2 * std::numeric_limits<T>::digits)), pol), -2 * std::numeric_limits<T>::digits);
    }

    expon = 1 + boost::math::ccmath::ilogb(val);
    
    if (T remain = boost::math::ccmath::scalbn(val, -expon); remain * std::numeric_limits<T>::radix == 1)
    {
        --expon; // when val is a power of two we must reduce the exponent
    }

    T diff = boost::math::ccmath::scalbn(static_cast<T>(1), expon - std::numeric_limits<T>::digits);
    if (diff == 0)
    {
        diff = detail::get_smallest_value<T>();
    }
    return val - diff;
} // float_prior_imp

} // namespace detail

template <typename T, typename Policy, typename result_type = tools::promote_args_t<T>>
constexpr result_type float_prior(const T& val, const Policy& pol)
{
    return detail::float_prior_imp(detail::normalize_value(static_cast<result_type>(val), typename detail::has_hidden_guard_digits<result_type>::type()), std::integral_constant<bool, !std::numeric_limits<result_type>::is_specialized || (std::numeric_limits<result_type>::radix == 2)>(), pol);
}

template <typename T, typename result_type = tools::promote_args_t<T>>
constexpr result_type float_prior(const T& val)
{
    return float_prior(val, policies::policy<>());
}

template <typename T, typename U, typename Policy, typename result_type = tools::promote_args_t<T, U>>
constexpr result_type nextafter(const T& val, const U& direction, const Policy& pol)
{
    if (BOOST_MATH_IS_CONSTANT_EVALUATED(val))
    {
        if (val < direction)
        {
            return boost::math::ccmath::float_next<result_type>(val, pol);
        }
        else if (val == direction)
        {
            // IEC 60559 recommends that from is returned whenever from == to. These functions return to instead, 
            // which makes the behavior around zero consistent: std::nextafter(-0.0, +0.0) returns +0.0 and 
            // std::nextafter(+0.0, -0.0) returns -0.0.
            return direction;
        }

        return boost::math::ccmath::float_prior<result_type>(val, pol);
    }
    else
    {
        using std::nextafter;
        return nextafter(static_cast<result_type>(val), static_cast<result_type>(direction));
    }
}

template <typename T, typename U, typename result_type = tools::promote_args_t<T, U>>
constexpr result_type nextafter(const T& val, const U& direction)
{
    return boost::math::ccmath::nextafter(val, direction, policies::policy<>());
}

constexpr float nextafterf(float val, float direction)
{
    return boost::math::ccmath::nextafter(val, direction);
}

constexpr long double nextafterl(long double val, long double direction)
{
    return boost::math::ccmath::nextafter(val, direction);
}

constexpr float nexttoward(float val, long double direction)
{
    return static_cast<float>(boost::math::ccmath::nextafter(val, direction));
}

constexpr float nexttowardf(float val, long double direction)
{
    return static_cast<float>(boost::math::ccmath::nextafter(val, direction));
}

constexpr double nexttoward(double val, long double direction)
{
    return static_cast<double>(boost::math::ccmath::nextafter(val, direction));
}

constexpr long double nexttoward(long double val, long double direction)
{
    return boost::math::ccmath::nextafter(val, direction);
}

constexpr long double nexttowardl(long double val, long double direction)
{
    return boost::math::ccmath::nextafter(val, direction);
}

template <typename T, typename result_type = tools::promote_args_t<T, long double>, typename return_type = std::conditional_t<std::is_integral_v<T>, double, result_type>>
constexpr return_type nexttoward(T val, long double direction)
{
    return static_cast<return_type>(boost::math::ccmath::nextafter(static_cast<result_type>(val), direction));
}

} // Namespaces

#endif // BOOST_MATH_SPECIAL_NEXT_HPP
