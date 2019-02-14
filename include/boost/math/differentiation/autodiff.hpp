//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DIFFERENTIATION_AUTODIFF_HPP
#define BOOST_MATH_DIFFERENTIATION_AUTODIFF_HPP

#include <boost/config.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/promotion.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <ostream>
#include <type_traits>

// Automatic Differentiation v1
namespace boost { namespace math { namespace differentiation { inline namespace autodiff_v1 {

namespace detail {

template<typename RealType, typename... RealTypes>
struct promote_args_n { using type = typename boost::math::tools::promote_args_2<RealType,
    typename promote_args_n<RealTypes...>::type>::type; };

template<typename RealType>
struct promote_args_n<RealType> { using type = typename boost::math::tools::promote_arg<RealType>::type; };

} // namespace detail

template<typename RealType, typename... RealTypes>
using promote = typename detail::promote_args_n<RealType,RealTypes...>::type;

namespace detail {

template<typename RealType, size_t Order>
class fvar;

template <typename>
struct get_depth : std::integral_constant<size_t, 0> {};

template <typename RealType, size_t Order>
struct get_depth<fvar<RealType,Order>> : std::integral_constant<size_t,get_depth<RealType>::value+1> {};

template <typename>
struct get_order_sum : std::integral_constant<size_t, 0> {};

template <typename RealType, size_t Order>
struct get_order_sum<fvar<RealType,Order>> : std::integral_constant<size_t,get_order_sum<RealType>::value+Order> {};

// Get non-fvar<> root type T of autodiff_fvar<T,O0,O1,O2,...>.
template<typename RealType>
struct get_root_type { using type = RealType; };

template<typename RealType, size_t Order>
struct get_root_type<fvar<RealType,Order>> { using type = typename get_root_type<RealType>::type; };

// Get type from descending Depth levels into fvar<>.
template<typename RealType, size_t Depth>
struct type_at { using type = RealType; };

template<typename RealType, size_t Order, size_t Depth>
struct type_at<fvar<RealType,Order>,Depth> { using type =
    typename std::conditional<Depth==0, fvar<RealType,Order>, typename type_at<RealType,Depth-1>::type>::type; };

template<typename RealType, size_t Depth>
using get_type_at = typename type_at<RealType,Depth>::type;

// Satisfies Boost's Conceptual Requirements for Real Number Types.
// https://www.boost.org/libs/math/doc/html/math_toolkit/real_concepts.html
template<typename RealType, size_t Order>
class fvar
{
    std::array<RealType,Order+1> v;

  public:

    using root_type = typename get_root_type<RealType>::type; // RealType in the root fvar<RealType,Order>.

    fvar() = default;

    // Initialize a variable or constant.
    fvar(const root_type&, const bool is_variable);

    // RealType(cr) | RealType | RealType is copy constructible.
    fvar(const fvar&) = default;

    // Be aware of implicit casting from one fvar<> type to another by this copy constructor.
    template<typename RealType2, size_t Order2>
    fvar(const fvar<RealType2,Order2>&);

    // RealType(ca) | RealType | RealType is copy constructible from the arithmetic types.
    explicit fvar(const root_type&); // Initialize a constant. (No epsilon terms.)

    template<typename RealType2>
    fvar(const RealType2& ca); // Supports any RealType2 for which static_cast<root_type>(ca) compiles.

    // r = cr | RealType& | Assignment operator.
    fvar& operator=(const fvar&) = default;

    // r = ca | RealType& | Assignment operator from the arithmetic types.
    // Handled by constructor that takes a single parameter of generic type.
    //fvar& operator=(const root_type&); // Set a constant.

    // r += cr | RealType& | Adds cr to r.
    template<typename RealType2, size_t Order2>
    fvar& operator+=(const fvar<RealType2,Order2>&);

    // r += ca | RealType& | Adds ar to r.
    fvar& operator+=(const root_type&);

    // r -= cr | RealType& | Subtracts cr from r.
    template<typename RealType2, size_t Order2>
    fvar& operator-=(const fvar<RealType2,Order2>&);

    // r -= ca | RealType& | Subtracts ca from r.
    fvar& operator-=(const root_type&);

    // r *= cr | RealType& | Multiplies r by cr.
    template<typename RealType2, size_t Order2>
    fvar& operator*=(const fvar<RealType2,Order2>&);

    // r *= ca | RealType& | Multiplies r by ca.
    fvar& operator*=(const root_type&);

    // r /= cr | RealType& | Divides r by cr.
    template<typename RealType2, size_t Order2>
    fvar& operator/=(const fvar<RealType2,Order2>&);

    // r /= ca | RealType& | Divides r by ca.
    fvar& operator/=(const root_type&);

    // -r | RealType | Unary Negation.
    fvar operator-() const;

    // +r | RealType& | Identity Operation.
    const fvar& operator+() const;

    // cr + cr2 | RealType | Binary Addition
    template<typename RealType2, size_t Order2>
    promote<fvar,fvar<RealType2,Order2>> operator+(const fvar<RealType2,Order2>&) const;

    // cr + ca | RealType | Binary Addition
    fvar operator+(const root_type&) const;

    // ca + cr | RealType | Binary Addition
    template<typename RealType2, size_t Order2>
    friend fvar<RealType2,Order2>
        operator+(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr - cr2 | RealType | Binary Subtraction
    template<typename RealType2, size_t Order2>
    promote<fvar,fvar<RealType2,Order2>> operator-(const fvar<RealType2,Order2>&) const;

    // cr - ca | RealType | Binary Subtraction
    fvar operator-(const root_type&) const;

    // ca - cr | RealType | Binary Subtraction
    template<typename RealType2, size_t Order2>
    friend fvar<RealType2,Order2>
        operator-(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr * cr2 | RealType | Binary Multiplication
    template<typename RealType2, size_t Order2>
    promote<fvar,fvar<RealType2,Order2>> operator*(const fvar<RealType2,Order2>&) const;

    // cr * ca | RealType | Binary Multiplication
    fvar operator*(const root_type&) const;

    // ca * cr | RealType | Binary Multiplication
    template<typename RealType2, size_t Order2>
    friend fvar<RealType2,Order2>
        operator*(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr / cr2 | RealType | Binary Subtraction
    template<typename RealType2, size_t Order2>
    promote<fvar,fvar<RealType2,Order2>> operator/(const fvar<RealType2,Order2>&) const;

    // cr / ca | RealType | Binary Subtraction
    fvar operator/(const root_type&) const;

    // ca / cr | RealType | Binary Subtraction
    template<typename RealType2, size_t Order2>
    friend fvar<RealType2,Order2>
        operator/(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr == cr2 | bool | Equality Comparison
    template<typename RealType2, size_t Order2> // This only compares the root term. All other terms are ignored.
    bool operator==(const fvar<RealType2,Order2>&) const;

    // cr == ca | bool | Equality Comparison
    bool operator==(const root_type&) const;

    // ca == cr | bool | Equality Comparison
    template<typename RealType2, size_t Order2> // This only compares the root term. All other terms are ignored.
    friend bool operator==(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr != cr2 | bool | Inequality Comparison
    template<typename RealType2, size_t Order2>
    bool operator!=(const fvar<RealType2,Order2>&) const;

    // cr != ca | bool | Inequality Comparison
    bool operator!=(const root_type&) const;

    // ca != cr | bool | Inequality Comparison
    template<typename RealType2, size_t Order2>
    friend bool operator!=(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr <= cr2 | bool | Less than equal to.
    template<typename RealType2, size_t Order2>
    bool operator<=(const fvar<RealType2,Order2>&) const;

    // cr <= ca | bool | Less than equal to.
    bool operator<=(const root_type&) const;

    // ca <= cr | bool | Less than equal to.
    template<typename RealType2, size_t Order2>
    friend bool operator<=(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr >= cr2 | bool | Greater than equal to.
    template<typename RealType2, size_t Order2>
    bool operator>=(const fvar<RealType2,Order2>&) const;

    // cr >= ca | bool | Greater than equal to.
    bool operator>=(const root_type&) const;

    // ca >= cr | bool | Greater than equal to.
    template<typename RealType2, size_t Order2>
    friend bool operator>=(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr < cr2 | bool | Less than comparison.
    template<typename RealType2, size_t Order2>
    bool operator<(const fvar<RealType2,Order2>&) const;

    // cr < ca | bool | Less than comparison.
    bool operator<(const root_type&) const;

    // ca < cr | bool | Less than comparison.
    template<typename RealType2, size_t Order2>
    friend bool operator<(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // cr > cr2 | bool | Greater than comparison.
    template<typename RealType2, size_t Order2>
    bool operator>(const fvar<RealType2,Order2>&) const;

    // cr > ca | bool | Greater than comparison.
    bool operator>(const root_type&) const;

    // ca > cr | bool | Greater than comparison.
    template<typename RealType2, size_t Order2>
    friend bool operator>(const typename fvar<RealType2,Order2>::root_type&, const fvar<RealType2,Order2>&);

    // Will throw std::out_of_range if Order < order.
    template<typename... Orders>
    get_type_at<RealType, sizeof...(Orders)> at(size_t order, Orders... orders) const;

    template<typename... Orders>
    get_type_at<fvar, sizeof...(Orders)> derivative(Orders... orders) const;

    fvar inverse() const; // Multiplicative inverse.

    static constexpr size_t depth = get_depth<fvar>::value; // Number of nested std::array<RealType,Order>.

    static constexpr size_t order_sum = get_order_sum<fvar>::value;

    explicit operator root_type() const; // Must be explicit, otherwise overloaded operators are ambiguous.

    fvar& set_root(const root_type&);

    // Use when function returns derivatives.
    fvar apply(const std::function<root_type(size_t)>&) const;

    // Use when function returns derivative(i)/factorial(i) (slightly more efficient than apply().)
    fvar apply_with_factorials(const std::function<root_type(size_t)>&) const;

    // Same as apply() but uses horner method. May be more accurate in some cases but not as good with inf derivatives.
    fvar apply_with_horner(const std::function<root_type(size_t)>&) const;

    // Same as apply_with_factorials() but uses horner method.
    fvar apply_with_horner_factorials(const std::function<root_type(size_t)>&) const;

private:

    RealType epsilon_inner_product(size_t z0, size_t isum0, size_t m0,
        const fvar& cr, size_t z1, size_t isum1, size_t m1, size_t j) const;

    fvar epsilon_multiply(size_t z0, size_t isum0, const fvar& cr, size_t z1, size_t isum1) const;

    fvar epsilon_multiply(size_t z0, size_t isum0, const root_type& ca) const;

    fvar inverse_apply() const;

    fvar& multiply_assign_by_root_type(bool is_root, const root_type&);

    template<typename RealType2, size_t Orders2>
    friend class fvar;

    template<typename RealType2, size_t Order2>
    friend std::ostream& operator<<(std::ostream&, const fvar<RealType2,Order2>&);

// C++11 Compatibility
#ifdef BOOST_NO_CXX17_IF_CONSTEXPR
    template<typename RootType>
    void fvar_cpp11(std::true_type, const RootType& ca, const bool is_variable);

    template<typename RootType>
    void fvar_cpp11(std::false_type, const RootType& ca, const bool is_variable);

    template<typename... Orders>
    get_type_at<RealType, sizeof...(Orders)> at_cpp11(std::true_type, size_t order, Orders... orders) const;

    template<typename... Orders>
    get_type_at<RealType, sizeof...(Orders)> at_cpp11(std::false_type, size_t order, Orders... orders) const;

    template<typename SizeType>
    fvar epsilon_multiply_cpp11(std::true_type,
        SizeType z0, size_t isum0, const fvar& cr, size_t z1, size_t isum1) const;

    template<typename SizeType>
    fvar epsilon_multiply_cpp11(std::false_type,
        SizeType z0, size_t isum0, const fvar& cr, size_t z1, size_t isum1) const;

    template<typename SizeType>
    fvar epsilon_multiply_cpp11(std::true_type, SizeType z0, size_t isum0, const root_type& ca) const;

    template<typename SizeType>
    fvar epsilon_multiply_cpp11(std::false_type, SizeType z0, size_t isum0, const root_type& ca) const;

    template<typename RootType>
    fvar& multiply_assign_by_root_type_cpp11(std::true_type, bool is_root, const RootType& ca);

    template<typename RootType>
    fvar& multiply_assign_by_root_type_cpp11(std::false_type, bool is_root, const RootType& ca);

    template<typename RootType>
    fvar& set_root_cpp11(std::true_type, const RootType& root);

    template<typename RootType>
    fvar& set_root_cpp11(std::false_type, const RootType& root);
#endif
};

// C++11 compatibility
#ifdef BOOST_NO_CXX17_IF_CONSTEXPR
#  define BOOST_AUTODIFF_IF_CONSTEXPR
#else
#  define BOOST_AUTODIFF_IF_CONSTEXPR constexpr
#endif

// Standard Library Support Requirements

// fabs(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> fabs(const fvar<RealType,Order>&);

// abs(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> abs(const fvar<RealType,Order>&);

// ceil(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> ceil(const fvar<RealType,Order>&);

// floor(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> floor(const fvar<RealType,Order>&);

// exp(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> exp(const fvar<RealType,Order>&);

// pow(cr, ca) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> pow(const fvar<RealType,Order>&,const typename fvar<RealType,Order>::root_type&);

// pow(ca, cr) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> pow(const typename fvar<RealType,Order>::root_type&,const fvar<RealType,Order>&);

// pow(cr1, cr2) | RealType
template<typename RealType1, size_t Order1, typename RealType2, size_t Order2>
promote<fvar<RealType1,Order1>,fvar<RealType2,Order2>>
    pow(const fvar<RealType1,Order1>&, const fvar<RealType2,Order2>&);

// sqrt(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> sqrt(const fvar<RealType,Order>&);

// log(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> log(const fvar<RealType,Order>&);

// frexp(cr1, &i) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> frexp(const fvar<RealType,Order>&, int*);

// ldexp(cr1, i) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> ldexp(const fvar<RealType,Order>&, int);

// cos(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> cos(const fvar<RealType,Order>&);

// sin(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> sin(const fvar<RealType,Order>&);

// asin(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> asin(const fvar<RealType,Order>&);

// tan(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> tan(const fvar<RealType,Order>&);

// atan(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> atan(const fvar<RealType,Order>&);

// fmod(cr1,cr2) | RealType
template<typename RealType1, size_t Order1, typename RealType2, size_t Order2>
promote<fvar<RealType1,Order1>,fvar<RealType2,Order2>>
    fmod(const fvar<RealType1,Order1>&, const fvar<RealType2,Order2>&);

// round(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> round(const fvar<RealType,Order>&);

// iround(cr1) | int
template<typename RealType, size_t Order>
int iround(const fvar<RealType,Order>&);

// trunc(cr1) | RealType
template<typename RealType, size_t Order>
fvar<RealType,Order> trunc(const fvar<RealType,Order>&);

// itrunc(cr1) | int
template<typename RealType, size_t Order>
int itrunc(const fvar<RealType,Order>&);

// Additional functions
template<typename RealType, size_t Order>
fvar<RealType,Order> acos(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> acosh(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> asinh(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> atanh(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> cosh(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> erf(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> erfc(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> lambert_w0(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> sinc(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> sinh(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
fvar<RealType,Order> tanh(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
long lround(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
long long llround(const fvar<RealType,Order>&);

template<typename RealType, size_t Order>
long double truncl(const fvar<RealType,Order>&);

// Compile-time test for fvar<> type.
template<typename>
struct is_fvar : std::false_type {};

template<typename RealType, size_t Order>
struct is_fvar<fvar<RealType,Order>> : std::true_type {};

template<typename RealType, size_t Order, size_t... Orders> // specialized for fvar<> below.
struct nest_fvar { using type = fvar<typename nest_fvar<RealType,Orders...>::type,Order>; };

template<typename RealType, size_t Order>
struct nest_fvar<RealType,Order> { using type = fvar<RealType,Order>; };

} // namespace detail

template<typename RealType, size_t Order, size_t... Orders>
using autodiff_fvar = typename detail::nest_fvar<RealType,Order,Orders...>::type;

template<typename RealType, size_t Order, size_t... Orders>
autodiff_fvar<RealType,Order,Orders...> make_fvar(const RealType& ca)
{
    return autodiff_fvar<RealType,Order,Orders...>(ca, true);
}

namespace detail {

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType, size_t Order>
fvar<RealType,Order>::fvar(const root_type& ca, const bool is_variable)
{
    if constexpr (is_fvar<RealType>::value)
    {
        v.front() = RealType(ca, is_variable);
        if constexpr (0 < Order)
            std::fill(v.begin()+1, v.end(), static_cast<RealType>(0));
    }
    else
    {
        v.front() = ca;
        if constexpr (0 < Order)
            v[1] = static_cast<root_type>(static_cast<int>(is_variable));
        if constexpr (1 < Order)
            std::fill(v.begin()+2, v.end(), static_cast<RealType>(0));
    }
}
#endif

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
fvar<RealType,Order>::fvar(const fvar<RealType2,Order2>& cr)
{
    for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
        v[i] = static_cast<RealType>(cr.v[i]);
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order2 < Order)
        std::fill(v.begin()+(Order2+1), v.end(), static_cast<RealType>(0));
}

template<typename RealType, size_t Order>
fvar<RealType,Order>::fvar(const root_type& ca)
:    v{{static_cast<RealType>(ca)}}
{
}

template<typename RealType, size_t Order>
template<typename RealType2>
fvar<RealType,Order>::fvar(const RealType2& ca)
:    v{{static_cast<RealType>(ca)}} // Can cause compiler error if RealType2 cannot be cast to root_type.
{
}

/*
template<typename RealType, size_t Order>
fvar<RealType,Order>& fvar<RealType,Order>::operator=(const root_type& ca)
{
    v.front() = static_cast<RealType>(ca);
    if BOOST_AUTODIFF_IF_CONSTEXPR (0 < Order)
        std::fill(v.begin()+1, v.end(), static_cast<RealType>(0));
    return *this;
}
*/

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
fvar<RealType,Order>& fvar<RealType,Order>::operator+=(const fvar<RealType2,Order2>& cr)
{
    for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
        v[i] += cr.v[i];
    return *this;
}

template<typename RealType, size_t Order>
fvar<RealType,Order>& fvar<RealType,Order>::operator+=(const root_type& ca)
{
    v.front() += ca;
    return *this;
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
fvar<RealType,Order>& fvar<RealType,Order>::operator-=(const fvar<RealType2,Order2>& cr)
{
    for (size_t i=0 ; i<=Order ; ++i)
        v[i] -= cr.v[i];
    return *this;
}

template<typename RealType, size_t Order>
fvar<RealType,Order>& fvar<RealType,Order>::operator-=(const root_type& ca)
{
    v.front() -= ca;
    return *this;
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
fvar<RealType,Order>& fvar<RealType,Order>::operator*=(const fvar<RealType2,Order2>& cr)
{
    const promote<RealType,RealType2> zero(0);
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order <= Order2)
        for (size_t i=0, j=Order ; i<=Order ; ++i, --j)
            v[j] = std::inner_product(v.cbegin(), v.cend()-i, cr.v.crbegin()+i, zero);
    else
    {
        for (size_t i=0, j=Order ; i<=Order-Order2 ; ++i, --j)
            v[j] = std::inner_product(cr.v.cbegin(), cr.v.cend(), v.crbegin()+i, zero);
        for (size_t i=Order-Order2+1, j=Order2-1 ; i<=Order ; ++i, --j)
            v[j] = std::inner_product(cr.v.cbegin(), cr.v.cbegin()+(j+1), v.crbegin()+i, zero);
    }
    return *this;
}

template<typename RealType, size_t Order>
fvar<RealType,Order>& fvar<RealType,Order>::operator*=(const root_type& ca)
{
    return multiply_assign_by_root_type(true, ca);
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
fvar<RealType,Order>& fvar<RealType,Order>::operator/=(const fvar<RealType2,Order2>& cr)
{
    const RealType zero(0);
    v.front() /= cr.v.front();
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
        for (size_t i=1, j=Order2-1, k=Order ; i<=Order ; ++i, --j, --k)
            (v[i] -= std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, v.crbegin()+k, zero)) /= cr.v.front();
    else if BOOST_AUTODIFF_IF_CONSTEXPR (0 < Order2)
        for (size_t i=1, j=Order2-1, k=Order ; i<=Order ; ++i, j&&--j, --k)
            (v[i] -= std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, v.crbegin()+k, zero)) /= cr.v.front();
    else
        for (size_t i=1 ; i<=Order ; ++i)
            v[i] /= cr.v.front();
    return *this;
}

template<typename RealType, size_t Order>
fvar<RealType,Order>& fvar<RealType,Order>::operator/=(const root_type& ca)
{
    std::for_each(v.begin(), v.end(), [&ca](RealType& x) { x /= ca; });
    return *this;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::operator-() const
{
    fvar<RealType,Order> retval;
    for (size_t i=0 ; i<=Order ; ++i)
        retval.v[i] = -v[i];
    return retval;
}

template<typename RealType, size_t Order>
const fvar<RealType,Order>& fvar<RealType,Order>::operator+() const
{
    return *this;
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
promote<fvar<RealType,Order>,fvar<RealType2,Order2>>
    fvar<RealType,Order>::operator+(const fvar<RealType2,Order2>& cr) const
{
    promote<fvar<RealType,Order>,fvar<RealType2,Order2>> retval;
    for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
        retval.v[i] = v[i] + cr.v[i];
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
        for (size_t i=Order+1 ; i<=Order2 ; ++i)
            retval.v[i] = cr.v[i];
    else if BOOST_AUTODIFF_IF_CONSTEXPR (Order2 < Order)
        for (size_t i=Order2+1 ; i<=Order ; ++i)
            retval.v[i] = v[i];
    return retval;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::operator+(const root_type& ca) const
{
    fvar<RealType,Order> retval(*this);
    retval.v.front() += ca;
    return retval;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> operator+(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return cr + ca;
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
promote<fvar<RealType,Order>,fvar<RealType2,Order2>>
    fvar<RealType,Order>::operator-(const fvar<RealType2,Order2>& cr) const
{
    promote<fvar<RealType,Order>,fvar<RealType2,Order2>> retval;
    for (size_t i=0 ; i<=std::min(Order,Order2) ; ++i)
        retval.v[i] = v[i] - cr.v[i];
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
        for (size_t i=Order+1 ; i<=Order2 ; ++i)
            retval.v[i] = -cr.v[i];
    else if BOOST_AUTODIFF_IF_CONSTEXPR (Order2 < Order)
        for (size_t i=Order2+1 ; i<=Order ; ++i)
            retval.v[i] = v[i];
    return retval;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::operator-(const root_type& ca) const
{
    fvar<RealType,Order> retval(*this);
    retval.v.front() -= ca;
    return retval;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> operator-(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return -cr += ca;
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
promote<fvar<RealType,Order>,fvar<RealType2,Order2>>
    fvar<RealType,Order>::operator*(const fvar<RealType2,Order2>& cr) const
{
    const promote<RealType,RealType2> zero(0);
    promote<fvar<RealType,Order>,fvar<RealType2,Order2>> retval;
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
        for (size_t i=0, j=Order, k=Order2 ; i<=Order2 ; ++i, j&&--j, --k)
            retval.v[i] = std::inner_product(v.cbegin(), v.cend()-j, cr.v.crbegin()+k, zero);
    else
        for (size_t i=0, j=Order2, k=Order ; i<=Order ; ++i, j&&--j, --k)
            retval.v[i] = std::inner_product(cr.v.cbegin(), cr.v.cend()-j, v.crbegin()+k, zero);
    return retval;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::operator*(const root_type& ca) const
{
    return fvar<RealType,Order>(*this) *= ca;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> operator*(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return cr * ca;
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
promote<fvar<RealType,Order>,fvar<RealType2,Order2>>
    fvar<RealType,Order>::operator/(const fvar<RealType2,Order2>& cr) const
{
    const promote<RealType,RealType2> zero(0);
    promote<fvar<RealType,Order>,fvar<RealType2,Order2>> retval;
    retval.v.front() = v.front() / cr.v.front();
    if BOOST_AUTODIFF_IF_CONSTEXPR (Order < Order2)
    {
        for (size_t i=1, j=Order2-1 ; i<=Order ; ++i, --j)
            retval.v[i] = (v[i] -
                std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, retval.v.crbegin()+(j+1), zero)) / cr.v.front();
        for (size_t i=Order+1, j=Order2-Order-1 ; i<=Order2 ; ++i, --j)
            retval.v[i] =
                -std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, retval.v.crbegin()+(j+1), zero) / cr.v.front();
    }
    else if BOOST_AUTODIFF_IF_CONSTEXPR (0 < Order2)
        for (size_t i=1, j=Order2-1, k=Order ; i<=Order ; ++i, j&&--j, --k)
            retval.v[i] =
                (v[i] - std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, retval.v.crbegin()+k, zero)) / cr.v.front();
    else
        for (size_t i=1 ; i<=Order ; ++i)
            retval.v[i] = v[i] / cr.v.front();
    return retval;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::operator/(const root_type& ca) const
{
    return fvar<RealType,Order>(*this) /= ca;
}

template<typename RealType, size_t Order>
fvar<RealType,Order> operator/(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    fvar<RealType,Order> retval;
    retval.v.front() = ca / cr.v.front();
    if BOOST_AUTODIFF_IF_CONSTEXPR (0 < Order)
    {
        const RealType zero(0);
        for (size_t i=1, j=Order-1 ; i<=Order ; ++i, --j)
            retval.v[i] = -std::inner_product(cr.v.cbegin()+1, cr.v.cend()-j, retval.v.crbegin()+(j+1), zero)
                / cr.v.front();
    }
    return retval;
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
bool fvar<RealType,Order>::operator==(const fvar<RealType2,Order2>& cr) const
{
    return v.front() == cr.v.front();
}

template<typename RealType, size_t Order>
bool fvar<RealType,Order>::operator==(const root_type& ca) const
{
    return v.front() == ca;
}

template<typename RealType, size_t Order>
bool operator==(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return ca == cr.v.front();
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
bool fvar<RealType,Order>::operator!=(const fvar<RealType2,Order2>& cr) const
{
    return v.front() != cr.v.front();
}

template<typename RealType, size_t Order>
bool fvar<RealType,Order>::operator!=(const root_type& ca) const
{
    return v.front() != ca;
}

template<typename RealType, size_t Order>
bool operator!=(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return ca != cr.v.front();
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
bool fvar<RealType,Order>::operator<=(const fvar<RealType2,Order2>& cr) const
{
    return v.front() <= cr.v.front();
}

template<typename RealType, size_t Order>
bool fvar<RealType,Order>::operator<=(const root_type& ca) const
{
    return v.front() <= ca;
}

template<typename RealType, size_t Order>
bool operator<=(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return ca <= cr.v.front();
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
bool fvar<RealType,Order>::operator>=(const fvar<RealType2,Order2>& cr) const
{
    return v.front() >= cr.v.front();
}

template<typename RealType, size_t Order>
bool fvar<RealType,Order>::operator>=(const root_type& ca) const
{
    return v.front() >= ca;
}

template<typename RealType, size_t Order>
bool operator>=(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return ca >= cr.v.front();
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
bool fvar<RealType,Order>::operator<(const fvar<RealType2,Order2>& cr) const
{
    return v.front() < cr.v.front();
}

template<typename RealType, size_t Order>
bool fvar<RealType,Order>::operator<(const root_type& ca) const
{
    return v.front() < ca;
}

template<typename RealType, size_t Order>
bool operator<(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return ca < cr.v.front();
}

template<typename RealType, size_t Order>
template<typename RealType2, size_t Order2>
bool fvar<RealType,Order>::operator>(const fvar<RealType2,Order2>& cr) const
{
    return v.front() > cr.v.front();
}

template<typename RealType, size_t Order>
bool fvar<RealType,Order>::operator>(const root_type& ca) const
{
    return v.front() > ca;
}

template<typename RealType, size_t Order>
bool operator>(const typename fvar<RealType,Order>::root_type& ca, const fvar<RealType,Order>& cr)
{
    return ca > cr.v.front();
}

/*** Other methods and functions ***/

// f : order -> derivative(order)
template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::apply(const std::function<root_type(size_t)>& f) const
{
    const fvar<RealType,Order> epsilon = fvar<RealType,Order>(*this).set_root(0);
    fvar<RealType,Order> epsilon_i = fvar<RealType,Order>(1); // epsilon to the power of i
    fvar<RealType,Order> accumulator = fvar<RealType,Order>(f(0));
    for (size_t i=1 ; i<=order_sum ; ++i)
    {    // accumulator += (epsilon_i *= epsilon) * (f(i) / boost::math::factorial<root_type>(i));
        epsilon_i = epsilon_i.epsilon_multiply(i-1, 0, epsilon, 1, 0);
        accumulator += epsilon_i.epsilon_multiply(i, 0, f(i) / boost::math::factorial<root_type>(i));
    }
    return accumulator;
}

// f : order -> derivative(order)/factorial(order)
// Use this when the computation of the derivatives already includes the factorial terms. E.g. See atan().
template<typename RealType, size_t Order>
fvar<RealType,Order>
    fvar<RealType,Order>::apply_with_factorials(const std::function<root_type(size_t)>& f) const
{
    const fvar<RealType,Order> epsilon = fvar<RealType,Order>(*this).set_root(0);
    fvar<RealType,Order> epsilon_i = fvar<RealType,Order>(1); // epsilon to the power of i
    fvar<RealType,Order> accumulator = fvar<RealType,Order>(f(0));
    for (size_t i=1 ; i<=order_sum ; ++i)
    {    // accumulator += (epsilon_i *= epsilon) * f(i);
        epsilon_i = epsilon_i.epsilon_multiply(i-1, 0, epsilon, 1, 0);
        accumulator += epsilon_i.epsilon_multiply(i, 0, f(i));
    }
    return accumulator;
}

// f : order -> derivative(order)
template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::apply_with_horner(const std::function<root_type(size_t)>& f) const
{
    const fvar<RealType,Order> epsilon = fvar<RealType,Order>(*this).set_root(0);
    fvar<RealType,Order> accumulator(static_cast<root_type>(f(order_sum)/boost::math::factorial<root_type>(order_sum)));
    for (size_t i=order_sum ; i-- ;)
        (accumulator *= epsilon) += f(i) / boost::math::factorial<root_type>(i);
    return accumulator;
}

// f : order -> derivative(order)/factorial(order)
// Use this when the computation of the derivatives already includes the factorial terms. E.g. See atan().
template<typename RealType, size_t Order>
fvar<RealType,Order>
    fvar<RealType,Order>::apply_with_horner_factorials(const std::function<root_type(size_t)>& f) const
{
    const fvar<RealType,Order> epsilon = fvar<RealType,Order>(*this).set_root(0);
    fvar<RealType,Order> accumulator(f(order_sum));
    for (size_t i=order_sum ; i-- ;)
        (accumulator *= epsilon) += f(i);
    return accumulator;
}

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
// Can throw "std::out_of_range: array::at: __n (which is 7) >= _Nm (which is 7)"
template<typename RealType, size_t Order>
template<typename... Orders>
get_type_at<RealType,sizeof...(Orders)> fvar<RealType,Order>::at(size_t order, Orders... orders) const
{
    if constexpr (0 < sizeof...(Orders))
        return v.at(order).at(orders...);
    else
        return v.at(order);
}
#endif

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
// Can throw "std::out_of_range: array::at: __n (which is 7) >= _Nm (which is 7)"
template<typename RealType, size_t Order>
template<typename... Orders>
get_type_at<fvar<RealType,Order>,sizeof...(Orders)> fvar<RealType,Order>::derivative(Orders... orders) const
{
    static_assert(sizeof...(Orders) <= depth, "Number of parameters to derivative(...) cannot exceed fvar::depth.");
    return at(orders...) * (... * boost::math::factorial<root_type>(orders));
}
#endif

template<typename RealType, size_t Order>
RealType fvar<RealType,Order>::epsilon_inner_product(size_t z0, size_t isum0, size_t m0,
    const fvar<RealType,Order>& cr, size_t z1, size_t isum1, size_t m1, size_t j) const
{
    static_assert(is_fvar<RealType>::value, "epsilon_inner_product() must have 1 < depth.");
    RealType accumulator = RealType();
    const size_t i0_max = m1 < j ? j-m1 : 0;
    for (size_t i0=m0, i1=j-m0 ; i0<=i0_max ; ++i0, --i1)
        accumulator += v.at(i0).epsilon_multiply(z0, isum0+i0, cr.v.at(i1), z1, isum1+i1);
    return accumulator;
}

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::epsilon_multiply(size_t z0, size_t isum0,
    const fvar<RealType,Order>& cr, size_t z1, size_t isum1) const
{
    const RealType zero(0);
    const size_t m0 = order_sum + isum0 < Order + z0 ? Order + z0 - (order_sum + isum0) : 0;
    const size_t m1 = order_sum + isum1 < Order + z1 ? Order + z1 - (order_sum + isum1) : 0;
    const size_t i_max = m0 + m1 < Order ? Order - (m0 + m1) : 0;
    fvar<RealType,Order> retval = fvar<RealType,Order>();
    if constexpr (is_fvar<RealType>::value)
        for (size_t i=0, j=Order ; i<=i_max ; ++i, --j)
            retval.v[j] = epsilon_inner_product(z0, isum0, m0, cr, z1, isum1, m1, j);
    else
        for (size_t i=0, j=Order ; i<=i_max ; ++i, --j)
            retval.v[j] = std::inner_product(v.cbegin()+m0, v.cend()-(i+m1), cr.v.crbegin()+(i+m0), zero);
    return retval;
}
#endif

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
// When called from outside this method, z0 should be non-zero. Otherwise if z0=0 then it will give an
// incorrect result of 0 when the root value is 0 and ca=inf, when instead the correct product is nan.
// If z0=0 then use the regular multiply operator*() instead.
template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::epsilon_multiply(size_t z0, size_t isum0,
    const root_type& ca) const
{
    fvar<RealType,Order> retval(*this);
    const size_t m0 = order_sum + isum0 < Order + z0 ? Order + z0 - (order_sum + isum0) : 0;
    if constexpr (is_fvar<RealType>::value)
        for (size_t i=m0 ; i<=Order ; ++i)
            retval.v[i] = retval.v[i].epsilon_multiply(z0, isum0+i, ca);
    else
        for (size_t i=m0 ; i<=Order ; ++i)
            if (retval.v[i] != static_cast<RealType>(0))
                retval.v[i] *= ca;
    return retval;
}
#endif

template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::inverse() const
{
    return operator root_type() == 0 ? inverse_apply() : 1 / *this;
}

// This gives log(0.0) = depth(1)(-inf,inf,-inf,inf,-inf,inf)
// 1 / *this: log(0.0) = depth(1)(-inf,inf,-inf,-nan,-nan,-nan)
template<typename RealType, size_t Order>
fvar<RealType,Order> fvar<RealType,Order>::inverse_apply() const
{
    root_type derivatives[order_sum+1]; // LCOV_EXCL_LINE This causes a false negative on lcov coverage test.
    const root_type x0 = static_cast<root_type>(*this);
    *derivatives = 1 / x0;
    for (size_t i=1 ; i<=order_sum ; ++i)
        derivatives[i] = -derivatives[i-1] * i / x0;
    return apply([&derivatives](size_t j) { return derivatives[j]; });
}

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType, size_t Order>
fvar<RealType,Order>& fvar<RealType,Order>::multiply_assign_by_root_type(bool is_root, const root_type& ca)
{
    auto itr = v.begin();
    if constexpr (is_fvar<RealType>::value)
    {
        itr->multiply_assign_by_root_type(is_root, ca);
        for (++itr ; itr!=v.end() ; ++itr)
            itr->multiply_assign_by_root_type(false, ca);
    }
    else
    {
        if (is_root || *itr != 0)
            *itr *= ca; // Skip multiplication of 0 by ca=inf to avoid nan. Exception: root value is always multiplied.
        for (++itr ; itr!=v.end() ; ++itr)
            if (*itr != 0)
                *itr *= ca;
    }
    return *this;
}
#endif

template<typename RealType, size_t Order>
fvar<RealType,Order>::operator root_type() const
{
    return static_cast<root_type>(v.front());
}

#ifndef BOOST_NO_CXX17_IF_CONSTEXPR
template<typename RealType, size_t Order>
fvar<RealType,Order>& fvar<RealType,Order>::set_root(const root_type& root)
{
    if constexpr (is_fvar<RealType>::value)
        v.front().set_root(root);
    else
        v.front() = root;
    return *this;
}
#endif

// Standard Library Support Requirements

template<typename RealType, size_t Order>
fvar<RealType,Order> fabs(const fvar<RealType,Order>& cr)
{
    const typename fvar<RealType,Order>::root_type zero(0);
    return cr < zero ? -cr
        : cr == zero ? fvar<RealType,Order>() // Canonical fabs'(0) = 0.
        : cr; // Propagate NaN.
}

template<typename RealType, size_t Order>
fvar<RealType,Order> abs(const fvar<RealType,Order>& cr)
{
    return fabs(cr);
}

template<typename RealType, size_t Order>
fvar<RealType,Order> ceil(const fvar<RealType,Order>& cr)
{
    using std::ceil;
    return fvar<RealType,Order>(ceil(static_cast<typename fvar<RealType,Order>::root_type>(cr)));
}

template<typename RealType, size_t Order>
fvar<RealType,Order> floor(const fvar<RealType,Order>& cr)
{
    using std::floor;
    return fvar<RealType,Order>(floor(static_cast<typename fvar<RealType,Order>::root_type>(cr)));
}

template<typename RealType, size_t Order>
fvar<RealType,Order> exp(const fvar<RealType,Order>& cr)
{
    using std::exp;
    using root_type = typename fvar<RealType,Order>::root_type;
    const root_type d0 = exp(static_cast<root_type>(cr));
    return cr.apply_with_horner([&d0](size_t) { return d0; });
}

template<typename RealType, size_t Order>
fvar<RealType,Order> pow(const fvar<RealType,Order>& x,const typename fvar<RealType,Order>::root_type& y)
{
    using std::pow;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    root_type derivatives[order+1];
    const root_type x0 = static_cast<root_type>(x);
    size_t i = 0;
    root_type coef = 1;
    for (; i<=order && coef!=0 ; ++i)
    {
        derivatives[i] = coef * pow(x0, y-i);
        coef *= y - i;
    }
    return x.apply([&derivatives,i](size_t j) { return j < i ? derivatives[j] : 0; });
}

template<typename RealType, size_t Order>
fvar<RealType,Order> pow(const typename fvar<RealType,Order>::root_type& x,const fvar<RealType,Order>& y)
{
    using std::log;
    return exp(y*log(x));
}

template<typename RealType1, size_t Order1, typename RealType2, size_t Order2>
promote<fvar<RealType1,Order1>,fvar<RealType2,Order2>>
    pow(const fvar<RealType1,Order1>& x, const fvar<RealType2,Order2>& y)
{
    return exp(y*log(x));
}

template<typename RealType, size_t Order>
fvar<RealType,Order> sqrt(const fvar<RealType,Order>& cr)
{
    using std::sqrt;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    root_type derivatives[order+1];
    const root_type x = static_cast<root_type>(cr);
    *derivatives = sqrt(x);
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(*derivatives);
    else
    {
        root_type numerator = 0.5;
        root_type powers = 1;
        derivatives[1] = numerator / *derivatives;
        for (size_t i=2 ; i<=order ; ++i)
        {
            numerator *= -0.5 * ((i<<1)-3);
            powers *= x;
            derivatives[i] = numerator / (powers * *derivatives);
        }
        return cr.apply([&derivatives](size_t i) { return derivatives[i]; });
    }
}

// Natural logarithm. If cr==0 then derivative(i) may have nans due to nans from inverse().
template<typename RealType, size_t Order>
fvar<RealType,Order> log(const fvar<RealType,Order>& cr)
{
    using std::log;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = log(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        const auto d1 = make_fvar<root_type,order-1>(static_cast<root_type>(cr)).inverse(); // log'(x) = 1 / x
        return cr.apply_with_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> frexp(const fvar<RealType,Order>& cr, int* exp)
{
    using std::frexp;
    using root_type = typename fvar<RealType,Order>::root_type;
    frexp(static_cast<root_type>(cr), exp);
    return cr * std::exp2(-*exp);
}

template<typename RealType, size_t Order>
fvar<RealType,Order> ldexp(const fvar<RealType,Order>& cr, int exp)
{
    return cr * std::exp2(exp);
}

template<typename RealType, size_t Order>
fvar<RealType,Order> cos(const fvar<RealType,Order>& cr)
{
    using std::cos;
    using std::sin;
    using root_type = typename fvar<RealType,Order>::root_type;
    const root_type d0 = cos(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (fvar<RealType,Order>::order_sum == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        const root_type d1 = -sin(static_cast<root_type>(cr));
        const root_type derivatives[4] { d0, d1, -d0, -d1 };
        return cr.apply_with_horner([&derivatives](size_t i) { return derivatives[i&3]; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> sin(const fvar<RealType,Order>& cr)
{
    using std::sin;
    using std::cos;
    using root_type = typename fvar<RealType,Order>::root_type;
    const root_type d0 = sin(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (fvar<RealType,Order>::order_sum == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        const root_type d1 = cos(static_cast<root_type>(cr));
        const root_type derivatives[4] { d0, d1, -d0, -d1 };
        return cr.apply_with_horner([&derivatives](size_t i) { return derivatives[i&3]; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> asin(const fvar<RealType,Order>& cr)
{
    using std::asin;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = asin(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        auto d1 = make_fvar<root_type,order-1>(static_cast<root_type>(cr)); // asin'(x) = 1 / sqrt(1-x*x).
        d1 = sqrt(1-(d1*=d1)).inverse(); // asin(1): d1 = depth(1)(inf,inf,-nan,-nan,-nan)
        //d1 = sqrt((1-(d1*=d1)).inverse()); // asin(1): d1 = depth(1)(inf,-nan,-nan,-nan,-nan)
        return cr.apply_with_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> tan(const fvar<RealType,Order>& cr)
{
    return sin(cr) / cos(cr);
}

template<typename RealType, size_t Order>
fvar<RealType,Order> atan(const fvar<RealType,Order>& cr)
{
    using std::atan;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = atan(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        auto d1 = make_fvar<root_type,order-1>(static_cast<root_type>(cr));
        d1 = ((d1*=d1)+=1).inverse(); // atan'(x) = 1 / (x*x+1).
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType1, size_t Order1, typename RealType2, size_t Order2>
promote<fvar<RealType1,Order1>,fvar<RealType2,Order2>>
    fmod(const fvar<RealType1,Order1>& cr1, const fvar<RealType2,Order2>& cr2)
{
    using std::trunc;
    const auto numer = static_cast<typename fvar<RealType1,Order1>::root_type>(cr1);
    const auto denom = static_cast<typename fvar<RealType2,Order2>::root_type>(cr2);
    return cr1 - cr2 * trunc(numer/denom);
}

template<typename RealType, size_t Order>
fvar<RealType,Order> round(const fvar<RealType,Order>& cr)
{
    using std::round;
    return fvar<RealType,Order>(round(static_cast<typename fvar<RealType,Order>::root_type>(cr)));
}

template<typename RealType, size_t Order>
int iround(const fvar<RealType,Order>& cr)
{
    using boost::math::iround;
    return iround(static_cast<typename fvar<RealType,Order>::root_type>(cr));
}

template<typename RealType, size_t Order>
fvar<RealType,Order> trunc(const fvar<RealType,Order>& cr)
{
    using std::trunc;
    return fvar<RealType,Order>(trunc(static_cast<typename fvar<RealType,Order>::root_type>(cr)));
}

template<typename RealType, size_t Order>
int itrunc(const fvar<RealType,Order>& cr)
{
    using boost::math::itrunc;
    return itrunc(static_cast<typename fvar<RealType,Order>::root_type>(cr));
}

template<typename RealType, size_t Order>
std::ostream& operator<<(std::ostream& out, const fvar<RealType,Order>& cr)
{
    out << "depth(" << cr.depth << ')';
    for (size_t i=0 ; i<cr.v.size() ; ++i)
        out << (i?',':'(') << cr.v[i];
    return out << ')';
}

// Additional functions

template<typename RealType, size_t Order>
fvar<RealType,Order> acos(const fvar<RealType,Order>& cr)
{
    using std::acos;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = acos(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        auto x = make_fvar<root_type,order-1>(static_cast<root_type>(cr));
        const auto d1 = -sqrt(1-(x*=x)).inverse(); // acos'(x) = -1 / sqrt(1-x*x).
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> acosh(const fvar<RealType,Order>& cr)
{
    using std::acosh;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = acosh(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        auto x = make_fvar<root_type,order-1>(static_cast<root_type>(cr));
        const auto d1 = sqrt((x*=x)-1).inverse(); // acosh'(x) = 1 / sqrt(x*x-1).
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> asinh(const fvar<RealType,Order>& cr)
{
    using std::asinh;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = asinh(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        auto x = make_fvar<root_type,order-1>(static_cast<root_type>(cr));
        const auto d1 = sqrt((x*=x)+1).inverse(); // asinh'(x) = 1 / sqrt(x*x+1).
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> atanh(const fvar<RealType,Order>& cr)
{
    using std::atanh;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = atanh(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        auto x = make_fvar<root_type,order-1>(static_cast<root_type>(cr));
        const auto d1 = (1-(x*=x)).inverse(); // atanh'(x) = 1 / (1-x*x)
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> cosh(const fvar<RealType,Order>& cr)
{
    using std::cosh;
    using std::sinh;
    using root_type = typename fvar<RealType,Order>::root_type;
    const root_type d0 = cosh(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (fvar<RealType,Order>::order_sum == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        const root_type derivatives[2] { d0, sinh(static_cast<root_type>(cr)) };
        return cr.apply_with_horner([&derivatives](size_t i) { return derivatives[i&1]; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> erf(const fvar<RealType,Order>& cr)
{
    using std::erf;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = erf(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        auto x = make_fvar<root_type,order-1>(static_cast<root_type>(cr));
        const auto d1 = 2*boost::math::constants::one_div_root_pi<root_type>()*exp(-(x*=x)); // 2/sqrt(pi)*exp(-x*x)
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> erfc(const fvar<RealType,Order>& cr)
{
    using std::erfc;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    const root_type d0 = erfc(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        auto x = make_fvar<root_type,order-1>(static_cast<root_type>(cr));
        const auto d1 = -2*boost::math::constants::one_div_root_pi<root_type>()*exp(-(x*=x)); // erfc'(x)=-erf'(x)
        return cr.apply_with_horner_factorials([&d0,&d1](size_t i) { return i ? d1.at(i-1)/i : d0; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> lambert_w0(const fvar<RealType,Order>& cr)
{
    using boost::math::lambert_w0;
    using std::exp;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    root_type derivatives[order+1];
    *derivatives = lambert_w0(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(*derivatives);
    else
    {
        const root_type expw = exp(*derivatives);
        derivatives[1] = 1 / (static_cast<root_type>(cr) + expw);
        if BOOST_AUTODIFF_IF_CONSTEXPR (order == 1)
            return cr.apply([&derivatives](size_t i) { return derivatives[i]; });
        else
        {
            root_type d1powers = derivatives[1] * derivatives[1];
            const root_type x = derivatives[1] * expw;
            derivatives[2] = d1powers * (-1 - x);
            std::array<root_type,order> coef {{ -1, -1 }}; // as in derivatives[2].
            for (size_t n=3 ; n<=order ; ++n)
            {
                coef[n-1] = coef[n-2] * -static_cast<root_type>(2*n-3);
                for (size_t j=n-2 ; j!=0 ; --j)
                    (coef[j] *= -static_cast<root_type>(n-1)) -= (n+j-2) * coef[j-1];
                coef[0] *= -static_cast<root_type>(n-1);
                d1powers *= derivatives[1];
                derivatives[n] = d1powers * std::accumulate(coef.crend()-(n-1), coef.crend(), coef[n-1],
                    [&x](const root_type& a, const root_type& b) { return a*x + b; });
            }
            return cr.apply([&derivatives](size_t i) { return derivatives[i]; });
        }
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> sinc(const fvar<RealType,Order>& cr)
{
    if (cr != 0)
        return sin(cr) / cr;
    using root_type = typename fvar<RealType,Order>::root_type;
    constexpr size_t order = fvar<RealType,Order>::order_sum;
    root_type taylor[order+1] { 1 }; // sinc(0) = 1
    if BOOST_AUTODIFF_IF_CONSTEXPR (order == 0)
        return fvar<RealType,Order>(*taylor);
    else
    {
        for (size_t n=2 ; n<=order ; n+=2)
            taylor[n] = (1-static_cast<int>(n&2)) / boost::math::factorial<root_type>(n+1);
        return cr.apply_with_factorials([&taylor](size_t i) { return taylor[i]; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> sinh(const fvar<RealType,Order>& cr)
{
    using std::sinh;
    using std::cosh;
    using root_type = typename fvar<RealType,Order>::root_type;
    const root_type d0 = sinh(static_cast<root_type>(cr));
    if BOOST_AUTODIFF_IF_CONSTEXPR (fvar<RealType,Order>::order_sum == 0)
        return fvar<RealType,Order>(d0);
    else
    {
        const root_type derivatives[2] { d0, cosh(static_cast<root_type>(cr)) };
        return cr.apply_with_horner([&derivatives](size_t i) { return derivatives[i&1]; });
    }
}

template<typename RealType, size_t Order>
fvar<RealType,Order> tanh(const fvar<RealType,Order>& cr)
{
    const fvar<RealType,Order> exp2cr = exp(cr*2);
    return (exp2cr - 1) /= (exp2cr + 1);
}

template<typename RealType, size_t Order>
long lround(const fvar<RealType,Order>& cr)
{
    using std::lround;
    return lround(static_cast<typename fvar<RealType,Order>::root_type>(cr));
}

template<typename RealType, size_t Order>
long long llround(const fvar<RealType,Order>& cr)
{
    using std::llround;
    return llround(static_cast<typename fvar<RealType,Order>::root_type>(cr));
}

template<typename RealType, size_t Order>
long double truncl(const fvar<RealType,Order>& cr)
{
    using std::truncl;
    return truncl(static_cast<typename fvar<RealType,Order>::root_type>(cr));
}

} } } } } // namespace boost::math::differentiation::autodiff_v1::detail

namespace std {

/// boost::math::tools::digits<RealType>() is handled by this std::numeric_limits<> specialization,
/// and similarly for max_value, min_value, log_max_value, log_min_value, and epsilon.
template <typename RealType, size_t Order>
class numeric_limits<boost::math::differentiation::detail::fvar<RealType,Order>>
    : public numeric_limits<typename boost::math::differentiation::detail::fvar<RealType,Order>::root_type>
{ };

} // namespace std

namespace boost { namespace math { namespace tools {

// See boost/math/tools/promotion.hpp
template <typename RealType0, size_t Order0, typename RealType1, size_t Order1>
struct promote_args_2<differentiation::detail::fvar<RealType0,Order0>,differentiation::detail::fvar<RealType1,Order1>>
{
    using type = differentiation::detail::fvar<typename promote_args_2<RealType0,RealType1>::type,
#ifndef BOOST_NO_CXX14_CONSTEXPR
        std::max(Order0,Order1)>;
#else
        Order0 < Order1 ? Order1 : Order0>;
#endif
};

template <typename RealType0, size_t Order0, typename RealType1>
struct promote_args_2<differentiation::detail::fvar<RealType0,Order0>,RealType1>
{
    using type = differentiation::detail::fvar<typename promote_args_2<RealType0,RealType1>::type,Order0>;
};

template <typename RealType0, typename RealType1, size_t Order1>
struct promote_args_2<RealType0,differentiation::detail::fvar<RealType1,Order1>>
{
    using type = differentiation::detail::fvar<typename promote_args_2<RealType0,RealType1>::type,Order1>;
};

} } } // namespace boost::math::tools

#ifdef BOOST_NO_CXX17_IF_CONSTEXPR
#include "autodiff_cpp11.hpp"
#endif

#endif // BOOST_MATH_DIFFERENTIATION_AUTODIFF_HPP
