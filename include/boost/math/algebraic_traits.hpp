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

/*                                      PURPOSE
 * 
 * A facility to provide convenient access to algebraic traits of types.
 * 
 * For example:
 * 
 *      algebraic_traits<std::multiplies, int>::zero_element() -> 0
 *      algebraic_traits<std::multiplies, int>::identity_element() -> 1
 *      algebraic_traits<std::plus, std::string>::identity_element() -> ""
 * 
 *      algebraic_traits<std::multiplies, int>::is_associative() -> true
 *      algebraic_traits<std::multiplies, double>::is_associative() -> false
 *      algebraic_traits<std::plus, std::string>::is_associative() -> true
 *      algebraic_traits<std::multiplies, double>::is_commutative() -> true
 *      algebraic_traits<std::plus, std::string>::is_commutative() -> false
 * 
 * Or alternatively:
 * 
 *      std::multiplies<int> op1;
 *      std::multiplies<double> op2;
 *      std::plus<std::string> op3;
 *      algebraic_traits::zero_element(op1) -> 0
 *      algebraic_traits::zero_element(op2) -> 0.0
 *      algebraic_traits::zero_element(op3) -> ""
 * 
 * 
 * Properties to consider:
 * 
 * Medial
 * Distributive
 * Commutative
 * Idempotent
 * Associative
 * Unital
 * Cancellative
 */

namespace boost {
namespace math {
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
    BOOST_STATIC_CONSTANT(bool, is_specialized = true);
    
    BOOST_STATIC_CONSTEXPR
    T zero_element() { return T(0); }
    
    BOOST_STATIC_CONSTEXPR
    T identity_element() { return T(1); }
};



// What makes a type associative in programming: exactness?


//  plus

// The inefficient way to multiply inherit?
template <typename T>
struct algebraic_traits<std::plus, T> : public algebraic_traits_commutativity<std::plus, T>, public algebraic_traits_associativity<std::plus, T>
{
    BOOST_STATIC_CONSTANT(bool, is_specialized = true);

    BOOST_STATIC_CONSTEXPR
    T identity_element() { return T(0); }
};


template <>
struct algebraic_traits<std::plus, std::string>
{
    BOOST_STATIC_CONSTANT(bool, is_specialized = true);

    // TODO: Can be constexpr? String literal in C++14?
    BOOST_STATIC_CONSTEXPR
    std::string identity_element() { return std::string(); }
    
    BOOST_STATIC_CONSTEXPR
    bool is_commutative() { return false; }

    BOOST_STATIC_CONSTEXPR
    bool is_associative() { return true; }
};


} // math
} // boost

#endif
