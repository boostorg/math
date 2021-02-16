//  Copyright Peter Dimov 2015.
//  Copyright Matt Borland 2021.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
//  Template metaprogramming classes and functions to replace MPL
//  Source: http://www.pdimov.com/cpp2/simple_cxx11_metaprogramming.html

#ifndef BOOST_MATH_TOOLS_MP
#define BOOST_MATH_TOOLS_MP

#include <type_traits>
#include <cstddef>

namespace boost { namespace math { namespace tools {

// Types:
// Typelist 
template<typename... T>
struct mp_list {};

// Size_t
template<std::size_t N> 
using mp_size_t = std::integral_constant<std::size_t, N>;

// Boolean
template<bool B>
using mp_bool = std::integral_constant<bool, B>;




namespace detail {

// Size
template<typename L> 
struct mp_size_impl {};

template<template<typename...> class L, typename... T> // Template template parameter must use class
struct mp_size_impl<L<T...>>
{
    using type = std::integral_constant<std::size_t, sizeof...(T)>;
};

// Front
template<typename L>
struct mp_front_impl {};

template<template<typename...> class L, typename T1, typename... T> 
struct mp_front_impl<L<T1, T...>>
{
    using type = T1;
};

// Back
template<typename L>
struct mp_back_impl {};

template<template<typename...> class L, typename... T, typename B>
struct mp_back_impl<L<T..., B>>
{
    using type = B;
};

// At
// TODO - Use tree based lookup for larger typelists
// http://odinthenerd.blogspot.com/2017/04/tree-based-lookup-why-kvasirmpl-is.html
template<typename L, std::size_t>
struct mp_at {};

template<template<typename...> class L, typename T0, typename... T>
struct mp_at<L<T0, T...>, 0>
{
    using type = T0;
};

template<template<typename...> class L, typename T0, typename T1, typename... T>
struct mp_at<L<T0, T1, T...>, 1>
{
    using type = T1;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename... T>
struct mp_at<L<T0, T1, T2, T...>, 2>
{
    using type = T2;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename... T>
struct mp_at<L<T0, T1, T2, T3, T...>, 3>
{
    using type = T3;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T...>, 4>
{
    using type = T4;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T5, T...>, 5>
{
    using type = T5;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
         typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T5, T6, T...>, 6>
{
    using type = T6;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
         typename T7, typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T5, T6, T7, T...>, 7>
{
    using type = T7;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
         typename T7, typename T8, typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T5, T6, T7, T8, T...>, 8>
{
    using type = T8;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
         typename T7, typename T8, typename T9, typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T...>, 9>
{
    using type = T9;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
         typename T7, typename T8, typename T9, typename T10, typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T...>, 10>
{
    using type = T10;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
         typename T7, typename T8, typename T9, typename T10, typename T11, typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T...>, 11>
{
    using type = T11;
};

template<template<typename...> class L, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6,
         typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename... T>
struct mp_at<L<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T...>, 12>
{
    using type = T12;
};

// Push back
template<typename L, typename... T> 
struct mp_push_back_impl {};

template<template<typename...> class L, typename... U, typename... T> 
struct mp_push_back_impl<L<U...>, T...>
{
    using type = L<U..., T...>;
};

// Push front
template<typename L, typename... T>
struct mp_push_front_impl {};

template<template<typename...> class L, typename... U, typename... T>
struct mp_push_front_impl<L<U...>, T...>
{
    using type = L<T..., U...>;
};

// Find if
template<typename Iter>
constexpr std::size_t find_index(const Iter first, const Iter last)
{
    return first == last || *first ? 0 : 1 + find_index(first + 1, last);
}

template<typename L, template<typename...> class P>
struct mp_find_if_impl {};

template<template<typename...> class L, template<typename...> class P> 
struct mp_find_if_impl<L<>, P>
{
    using type = mp_size_t<0>;
};

template<template<typename...> class L, typename... T, template<typename...> class P> 
struct mp_find_if_impl<L<T...>, P>
{
    static constexpr bool v[] = {P<T>::value...};
    using type = mp_size_t<find_index(v, v + sizeof...(T))>;
};

} // namespace detail




// User Interfaces
template<typename T> 
using mp_size = typename detail::mp_size_impl<T>::type;

template<typename T>
using mp_front = typename detail::mp_front_impl<T>::type;

template<typename T>
using mp_back = typename detail::mp_back_impl<T>::type;

template<typename L, typename... T>
using mp_push_back = typename detail::mp_push_back_impl<L, T...>::type;

template<typename L, typename... T> 
using mp_push_front = typename detail::mp_push_front_impl<L, T...>::type;

template<typename L, template<typename...> class P> 
using mp_find_if = typename detail::mp_find_if_impl<L, P>::type;

}}} // namespaces

#endif // BOOST_MATH_TOOLS_MP
