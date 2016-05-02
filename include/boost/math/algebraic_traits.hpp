//  (C) Copyright Jeremy William Murphy 2016.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_ALGEBRAIC_TRAITS
#define BOOST_MATH_ALGEBRAIC_TRAITS

#include <boost/config.hpp>
#include <boost/core/enable_if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/type_traits/is_floating_point.hpp>

#include <limits>
#include <functional>
#include <string>

namespace boost {
namespace detail {
    
struct algebraic_traits_base
{
    BOOST_STATIC_CONSTANT(bool, is_specialized = false);
};

} // detail


template <template <typename> class BinaryOperation, typename T, typename Enable = void>
struct algebraic_traits : public detail::algebraic_traits_base
{
    // To provide default implementations or not?
};


// associativity

template <template <typename> class BinaryOperation, typename T, typename Enable = void>
struct algebraic_traits_associativity
{};


template <template <typename> class BinaryOperation, typename T>
struct algebraic_traits_associativity<BinaryOperation, T, 
typename enable_if< typename mpl::or_<is_integral<T>, mpl::bool_<std::numeric_limits<T>::is_exact> >::type >::type >
{
    BOOST_STATIC_CONSTEXPR
    bool is_associative() { return true; }
};


template <template <typename> class BinaryOperation, typename T>
struct algebraic_traits_associativity<BinaryOperation, T, 
typename enable_if< typename mpl::or_<is_floating_point<T>, mpl::bool_<!std::numeric_limits<T>::is_exact> >::type >::type >
{
    BOOST_STATIC_CONSTEXPR
    bool is_associative() { return false; }
};


// commutativity

template <template <typename> class BinaryOperation, typename T, typename Enable = void>
struct algebraic_traits : public detail::algebraic_traits_base
{
    // To provide default implementations or not?
};


// associativity

template <template <typename> class BinaryOperation, typename T, typename Enable = void>
struct algebraic_traits_commutativity
{};


template <typename T>
struct algebraic_traits_commutativity<std::plus, T>
{
    BOOST_STATIC_CONSTEXPR
    bool is_commutative() { return true; }
};


template <>
struct algebraic_traits_commutativity<std::plus, std::string>
{
    BOOST_STATIC_CONSTEXPR
    bool is_commutative() { return false; }
};


//  multiplication

template <typename T>
struct algebraic_traits< std::multiplies, T >
{
    BOOST_STATIC_CONSTEXPR
    T zero_element() { return T(0); }
    
    BOOST_STATIC_CONSTEXPR
    T identity_element() { return T(1); }
};



// What makes a type associative in programming: exactness?


//  plus

template <typename T>
struct algebraic_traits<std::plus, T> : public algebraic_traits_commutativity<std::plus, T>, public algebraic_traits_associativity<std::plus, T>
{
    BOOST_STATIC_CONSTEXPR
    T identity_element() { return T(0); }
};


template <>
struct algebraic_traits<std::plus, std::string>
{
    // TODO: Can be constexpr?
    BOOST_STATIC_CONSTEXPR
    std::string identity_element() { return std::string(); }
    
    BOOST_STATIC_CONSTEXPR
    bool is_commutative() { return false; }

    BOOST_STATIC_CONSTEXPR
    bool is_associative() { return true; }
};



} // boost

#endif
