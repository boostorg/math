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

} // namespace detail



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

}}} // namespaces

#endif // BOOST_MATH_TOOLS_MP
