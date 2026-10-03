//  (C) Copyright John Maddock 2010.
//  (C) Copyright Matt Borland 2024.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TUPLE_HPP_INCLUDED
#define BOOST_MATH_TUPLE_HPP_INCLUDED

#include <boost/math/tools/config.hpp>

#ifdef BOOST_MATH_ENABLE_CUDA

#include <boost/math/tools/type_traits.hpp>
#include <cuda/std/utility>
#include <cuda/std/tuple>

BOOST_MATH_NAMESPACE_BEGIN

BOOST_MATH_EXPORT using cuda::std::pair;
BOOST_MATH_EXPORT using cuda::std::tuple;

BOOST_MATH_EXPORT using cuda::std::make_pair;

BOOST_MATH_EXPORT using cuda::std::tie;
BOOST_MATH_EXPORT using cuda::std::get;

BOOST_MATH_EXPORT using cuda::std::tuple_size;
BOOST_MATH_EXPORT using cuda::std::tuple_element;

namespace detail {

template <typename T>
BOOST_MATH_GPU_ENABLED T&& forward(BOOST_MATH_NAMESPACE::remove_reference_t<T>& arg) noexcept
{
    return static_cast<T&&>(arg);
}

template <typename T>
BOOST_MATH_GPU_ENABLED T&& forward(BOOST_MATH_NAMESPACE::remove_reference_t<T>&& arg) noexcept
{
    static_assert(!BOOST_MATH_NAMESPACE::is_lvalue_reference<T>::value, "Cannot forward an rvalue as an lvalue.");
    return static_cast<T&&>(arg);
}

} // namespace detail

template <typename T, typename... Ts>
BOOST_MATH_GPU_ENABLED auto make_tuple(T&& t, Ts&&... ts) 
{
    return cuda::std::tuple<BOOST_MATH_NAMESPACE::decay_t<T>, BOOST_MATH_NAMESPACE::decay_t<Ts>...>(
        BOOST_MATH_NAMESPACE::detail::forward<T>(t), BOOST_MATH_NAMESPACE::detail::forward<Ts>(ts)...
    );
}

BOOST_MATH_NAMESPACE_END

#else

#ifndef BOOST_MATH_BUILD_MODULE
#include <tuple>
#endif

BOOST_MATH_NAMESPACE_BEGIN

BOOST_MATH_EXPORT using ::std::tuple;
BOOST_MATH_EXPORT using ::std::pair;

// [6.1.3.2] Tuple creation functions
BOOST_MATH_EXPORT using ::std::ignore;
BOOST_MATH_EXPORT using ::std::make_tuple;
BOOST_MATH_EXPORT using ::std::tie;
BOOST_MATH_EXPORT using ::std::get;

// [6.1.3.3] Tuple helper classes
BOOST_MATH_EXPORT using ::std::tuple_size;
BOOST_MATH_EXPORT using ::std::tuple_element;

// Pair helpers
BOOST_MATH_EXPORT using ::std::make_pair;

BOOST_MATH_NAMESPACE_END

#endif // BOOST_MATH_ENABLE_CUDA

#endif
