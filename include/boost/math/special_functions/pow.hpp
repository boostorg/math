//   Boost pow.hpp header file
//   Computes a power with exponent known at compile-time

//  (C) Copyright Bruno Lalande 2008.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org for updates, documentation, and revision history.


#ifndef BOOST_MATH_POW_HPP
#define BOOST_MATH_POW_HPP


#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/math/policies/policy.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/tools/promotion.hpp>
#include <boost/mpl/greater_equal.hpp>


namespace boost {
namespace math {

#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable:4702) // Unreachable code, only triggered in release mode and /W4
#endif

namespace detail {

   template <class T>
   inline BOOST_GPU_ENABLED BOOST_CONSTEXPR T square(const T& val) { return val * val; }


template <int N, int M = N%2>
struct positive_power
{
    template <typename T>
    static BOOST_GPU_ENABLED BOOST_CXX14_CONSTEXPR T result(T base)
    {
        return square(positive_power<N / 2>::result(base));
    }
};

template <int N>
struct positive_power<N, 1>
{
    template <typename T>
    static BOOST_GPU_ENABLED BOOST_CONSTEXPR T result(T base)
    {
        return base * square(positive_power<N / 2>::result(base));
    }
};

template <>
struct positive_power<1, 1>
{
    template <typename T>
    static BOOST_GPU_ENABLED BOOST_CONSTEXPR T result(T base){ return base; }
};


template <int N, bool>
struct power_if_positive
{
    template <typename T, class Policy>
    static BOOST_GPU_ENABLED BOOST_CONSTEXPR T result(T base, const Policy&)
    { return positive_power<N>::result(base); }
};

template <int N>
struct power_if_positive<N, false>
{
    template <typename T, class Policy>
    static BOOST_GPU_ENABLED BOOST_CONSTEXPR T result(T base, const Policy& policy)
    {
        return base == 0 ? policies::raise_overflow_error<T>(
                       "boost::math::pow(%1%)",
                       "Attempted to compute a negative power of 0",
                       policy
                   ) : T(1) / positive_power<-N>::result(base);
    }
};

template <>
struct power_if_positive<0, true>
{
    template <typename T, class Policy>
    static BOOST_GPU_ENABLED BOOST_CONSTEXPR T result(T base, const Policy& policy)
    {
        return base == 0 ? policies::raise_indeterminate_result_error<T>(
                       "boost::math::pow(%1%)",
                       "The result of pow<0>(%1%) is undetermined",
                       base,
                       T(1),
                       policy
                   ) : T(1);
    }
};


template <int N>
struct select_power_if_positive
{
    typedef typename mpl::greater_equal<
                         mpl::int_<N>,
                         mpl::int_<0>
                     >::type is_positive;

    typedef power_if_positive<N, is_positive::value> type;
};


}  // namespace detail


template <int N, typename T, class Policy>
inline BOOST_GPU_ENABLED BOOST_CONSTEXPR typename tools::promote_args<T>::type pow(T base, const Policy& policy)
{ 
   return detail::select_power_if_positive<N>::type::result(static_cast<typename tools::promote_args<T>::type>(base), policy);
}


template <int N, typename T>
inline BOOST_GPU_ENABLED BOOST_CONSTEXPR typename tools::promote_args<T>::type pow(T base)
{ return pow<N>(base, policies::policy<>()); }

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

}  // namespace math
}  // namespace boost


#endif
