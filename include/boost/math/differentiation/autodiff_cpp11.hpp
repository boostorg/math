//           Copyright Matthew Pulver 2018 - 2019.
// Distributed under the Boost Software License, Version 1.0.
//      (See accompanying file LICENSE_1_0.txt or copy at
//           https://www.boost.org/LICENSE_1_0.txt)

// Contributors:
//  * Kedar R. Bhat - C++11 compatibility.

// Notes:
//  * Any changes to this file should always be downstream from autodiff.cpp.
//    C++17 is a higher-level language and is easier to maintain. For example, a number of functions which are
//    lucidly read in autodiff.cpp are forced to be split into multiple structs/functions in this file for C++11.
//  * Use of typename RootType and SizeType is a hack to prevent Visual Studio 2015 from compiling functions
//    that are never called, that would otherwise produce compiler errors.

#ifndef BOOST_MATH_AUTODIFF_HPP
#   error "Do not #include this file directly. This should only be #included by autodiff.hpp for C++11 compatibility."
#endif

namespace boost { namespace math { namespace differentiation { namespace autodiff { inline namespace v1 {

template<typename RealType,size_t Order>
template<typename... Orders>
typename type_at<RealType, sizeof...(Orders)>::type
    dimension<RealType,Order>::at_cpp11(std::true_type, size_t order, Orders... orders) const
{
    return v.at(order);
}

template<typename RealType,size_t Order>
template<typename... Orders>
typename type_at<RealType, sizeof...(Orders)>::type
    dimension<RealType,Order>::at_cpp11(std::false_type, size_t order, Orders... orders) const
{
    return v.at(order).at(orders...);
}

// Can throw "std::out_of_range: array::at: __n (which is 7) >= _Nm (which is 7)"
template<typename RealType,size_t Order>
template<typename... Orders>
typename type_at<RealType,sizeof...(Orders)>::type dimension<RealType,Order>::at(size_t order, Orders... orders) const
{
    return at_cpp11(std::integral_constant<bool,sizeof...(orders)==0>{}, order, orders...);
}

template <typename>
struct depth_cpp11 : std::integral_constant<size_t, 0> {};

template <typename RealType,size_t Order>
struct depth_cpp11<dimension<RealType,Order>> : std::integral_constant<size_t,depth_cpp11<RealType>::value+1> {};

template<typename RealType,size_t Order>
constexpr size_t dimension<RealType,Order>::depth()
{
    return depth_cpp11<dimension<RealType,Order>>::value;
}

template <typename>
struct order_sum_cpp11 : std::integral_constant<size_t, 0> {};

template <typename RealType,size_t Order>
struct order_sum_cpp11<dimension<RealType,Order>> :
    std::integral_constant<size_t,order_sum_cpp11<RealType>::value+Order> {};

template<typename RealType,size_t Order>
constexpr size_t dimension<RealType,Order>::order_sum()
{
    return order_sum_cpp11<dimension<RealType,Order>>::value;
}

template<typename T, typename... Ts>
constexpr T product(Ts... factors)
{
    return 1;
}

template<typename T, typename... Ts>
constexpr T product(T factor, Ts... factors)
{
    return factor * product<T>(factors...);
}

// Can throw "std::out_of_range: array::at: __n (which is 7) >= _Nm (which is 7)"
template<typename RealType,size_t Order>
template<typename... Orders>
typename type_at<RealType,sizeof...(Orders)-1>::type dimension<RealType,Order>::derivative(Orders... orders) const
{
    static_assert(sizeof...(Orders) <= depth(),
        "Number of parameters to derivative(...) cannot exceed the number of dimensions in the dimension<...>.");
    return at(orders...) * product(boost::math::factorial<root_type>(orders)...);
}

template<typename RealType,size_t Order>
template<typename SizeType>
dimension<RealType,Order> dimension<RealType,Order>::epsilon_multiply_cpp11(std::true_type,
    SizeType z0, size_t isum0, const dimension<RealType,Order>& cr, size_t z1, size_t isum1) const
{
    const RealType zero{0};
    const size_t m0 = order_sum() + isum0 < Order + z0 ? Order + z0 - (order_sum() + isum0) : 0;
    const size_t m1 = order_sum() + isum1 < Order + z1 ? Order + z1 - (order_sum() + isum1) : 0;
    const size_t i_max = m0 + m1 < Order ? Order - (m0 + m1) : 0;
    dimension<RealType,Order> retval = dimension<RealType,Order>();
    for (size_t i=0, j=Order ; i<=i_max ; ++i, --j)
        retval.v[j] = epsilon_inner_product(z0, isum0, m0, cr, z1, isum1, m1, j);
    return retval;
}

template<typename RealType,size_t Order>
template<typename SizeType>
dimension<RealType,Order> dimension<RealType,Order>::epsilon_multiply_cpp11(std::false_type,
    SizeType z0, size_t isum0, const dimension<RealType,Order>& cr, size_t z1, size_t isum1) const
{
    const RealType zero{0};
    const size_t m0 = order_sum() + isum0 < Order + z0 ? Order + z0 - (order_sum() + isum0) : 0;
    const size_t m1 = order_sum() + isum1 < Order + z1 ? Order + z1 - (order_sum() + isum1) : 0;
    const size_t i_max = m0 + m1 < Order ? Order - (m0 + m1) : 0;
    dimension<RealType,Order> retval = dimension<RealType,Order>();
    for (size_t i=0, j=Order ; i<=i_max ; ++i, --j)
        retval.v[j] = std::inner_product(v.cbegin()+m0, v.cend()-(i+m1), cr.v.crbegin()+(i+m0), zero);
    return retval;
}

template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::epsilon_multiply(size_t z0, size_t isum0,
    const dimension<RealType,Order>& cr, size_t z1, size_t isum1) const
{
    return epsilon_multiply_cpp11(std::integral_constant<bool,is_dimension<RealType>::value>{},
        z0, isum0, cr, z1, isum1);
}

template<typename RealType,size_t Order>
template<typename SizeType>
dimension<RealType,Order> dimension<RealType,Order>::epsilon_multiply_cpp11(std::true_type,
    SizeType z0, size_t isum0, const root_type& ca) const
{
    dimension<RealType,Order> retval(*this);
    const size_t m0 = order_sum() + isum0 < Order + z0 ? Order + z0 - (order_sum() + isum0) : 0;
    for (size_t i=m0 ; i<=Order ; ++i)
        retval.v[i] = retval.v[i].epsilon_multiply(z0, isum0+i, ca);
    return retval;
}

template<typename RealType,size_t Order>
template<typename SizeType>
dimension<RealType,Order> dimension<RealType,Order>::epsilon_multiply_cpp11(std::false_type,
    SizeType z0, size_t isum0, const root_type& ca) const
{
    dimension<RealType,Order> retval(*this);
    const size_t m0 = order_sum() + isum0 < Order + z0 ? Order + z0 - (order_sum() + isum0) : 0;
    for (size_t i=m0 ; i<=Order ; ++i)
        if (retval.v[i] != static_cast<RealType>(0))
            retval.v[i] *= ca;
    return retval;
}

template<typename RealType,size_t Order>
dimension<RealType,Order> dimension<RealType,Order>::epsilon_multiply(size_t z0, size_t isum0,
    const root_type& ca) const
{
    return epsilon_multiply_cpp11(std::integral_constant<bool,is_dimension<RealType>::value>{}, z0, isum0, ca);
}

template<typename RealType,size_t Order>
template<typename RootType>
dimension<RealType,Order>& dimension<RealType,Order>::multiply_assign_by_root_type_cpp11(std::true_type,
    bool is_root, const RootType& ca)
{
    auto itr = v.begin();
    itr->multiply_assign_by_root_type(is_root, ca);
    for (++itr ; itr!=v.end() ; ++itr)
        itr->multiply_assign_by_root_type(false, ca);
    return *this;
}

template<typename RealType,size_t Order>
template<typename RootType>
dimension<RealType,Order>& dimension<RealType,Order>::multiply_assign_by_root_type_cpp11(std::false_type,
    bool is_root, const RootType& ca)
{
    auto itr = v.begin();
    if (is_root || *itr != 0)
        *itr *= ca; // Skip multiplication of 0 by ca=inf to avoid nan. Exception: root value is always multiplied.
    for (++itr ; itr!=v.end() ; ++itr)
        if (*itr != 0)
            *itr *= ca;
    return *this;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::multiply_assign_by_root_type(bool is_root, const root_type& ca)
{
    return multiply_assign_by_root_type_cpp11(std::integral_constant<bool,is_dimension<RealType>::value>{},
        is_root, ca);
}

template<typename RealType,size_t Order>
template<typename RootType>
dimension<RealType,Order>& dimension<RealType,Order>::set_root_cpp11(std::true_type, const RootType& root)
{
    v.front().set_root(root);
    return *this;
}

template<typename RealType,size_t Order>
template<typename RootType>
dimension<RealType,Order>& dimension<RealType,Order>::set_root_cpp11(std::false_type, const RootType& root)
{
    v.front() = root;
    return *this;
}

template<typename RealType,size_t Order>
dimension<RealType,Order>& dimension<RealType,Order>::set_root(const root_type& root)
{
    return set_root_cpp11(std::integral_constant<bool,is_dimension<RealType>::value>{}, root);
}

} } } } } // namespace boost::math::differentiation::autodiff::v1
