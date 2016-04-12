//  (C) Copyright Jeremy William Murphy 2016.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_COMMON_FACTOR_RT_HPP
#define BOOST_MATH_COMMON_FACTOR_RT_HPP

#include <boost/assert.hpp>
#include <boost/core/enable_if.hpp>
#include <boost/mpl/and.hpp>
#include <boost/type_traits.hpp>

#include <algorithm>
#include <limits>

namespace boost {
namespace math{

namespace detail
{
    template <typename T>
    inline bool odd(T const &x)
    {
        return static_cast<bool>(x & 0x1u);
    }
    
    
    template <typename T>
    inline bool even(T const &x)
    {
        return !odd(x);
    }
    
    
    /** Stein gcd (aka 'binary gcd')
     * 
     * From Mathematics to Generic Programming, Alexander Stepanov, Daniel Rose
     */
    template <typename SteinDomain>
    SteinDomain Stein_gcd(SteinDomain m, SteinDomain n)
    {
        using std::swap;
        BOOST_ASSERT(m >= 0);
        BOOST_ASSERT(n >= 0);
        if (m == SteinDomain(0))
            return n;
        if (n == SteinDomain(0))
            return m;
        // m > 0 && n > 0
        int d_m = 0;
        while (even(m))
        {
            m >>= 1;
            ++d_m;
        }
        int d_n = 0;
        while (even(n))
        {
            n >>= 1;
            ++d_n;
        }
        // odd(m) && odd(n)
        while (m != n)
        {
            if (n > m)
                swap(n, m);
            m -= n;
            do
                m >>= 1;
            while (even(m));
        }
        // m == n
        m <<= std::min(d_m, d_n);
        return m;
    }

    
    /** Euclidean algorithm
     * 
     * From Mathematics to Generic Programming, Alexander Stepanov, Daniel Rose
     * 
     */
    template <typename EuclideanDomain>
    inline EuclideanDomain Euclid_gcd(EuclideanDomain a, EuclideanDomain b)
    {
        using std::swap;
        while (b != EuclideanDomain(0))
        {
            a %= b;
            swap(a, b);
        }
        return a;
    }


    template <typename T>
    inline BOOST_DEDUCED_TYPENAME enable_if_c<mpl::and_< mpl::and_< has_right_shift_assign<T>, has_left_shift_assign<T> >, has_less<T> >::value, T>::type
       optimal_gcd_select(T const &a, T const &b)
    {
       return detail::Stein_gcd(a, b);
    }


    template <typename T>
    inline BOOST_DEDUCED_TYPENAME disable_if_c<mpl::and_< mpl::and_< has_right_shift_assign<T>, has_left_shift_assign<T> >, has_less<T> >::value, T>::type
       optimal_gcd_select(T const &a, T const &b)
    {
       return detail::Euclid_gcd(a, b);
    }

    //
    // To figure out whether a type has abs support or not, we use tribool logic:
    // 1) Definitely has abs,
    // 2) Definitely does not have abs,
    // 3) Don't know.
    // We want to call (3) only when all else fails as it doesn't work for built in types :(
    //
    typedef char nt[29761];
    typedef nt& no_type;
    template <class T> no_type abs(T const& ...);
    using std::abs;
    //
    // default template - don't know
    //
    template <class T, bool b, bool c>
    struct has_abs_imp
    {
       static const T val;
       BOOST_STATIC_CONSTANT(bool, value = sizeof(abs(val)) != sizeof(no_type));
    };
    //
    // Definitely no abs:
    //
    template <class T>
    struct has_abs_imp<T, true, false>
    {
       BOOST_STATIC_CONSTANT(bool, value = false);
    };
    //
    // Definitely has abs:
    //
    template <class T>
    struct has_abs_imp<T, false, true>
    {
       BOOST_STATIC_CONSTANT(bool, value = true);
    };

    template <class T>
    struct has_abs
    {
       BOOST_STATIC_CONSTANT(bool, value = (has_abs_imp<T, (boost::is_unsigned<T>::value || (std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_signed)), (boost::is_signed<T>::value || (std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::is_signed))>::value));
    };

} // namespace detail



template <typename Integer>
inline BOOST_DEDUCED_TYPENAME enable_if_c<detail::has_abs<Integer>::value, Integer>::type
gcd(Integer const &a, Integer const &b)
{
    using std::abs;
    return detail::optimal_gcd_select(static_cast<Integer>(abs(a)), static_cast<Integer>(abs(b)));
}


template <typename Integer>
inline BOOST_DEDUCED_TYPENAME disable_if_c<detail::has_abs<Integer>::value, Integer>::type
gcd(Integer const &a, Integer const &b)
{
    return detail::optimal_gcd_select(a, b);
}

namespace detail
{
   template <class T>
   inline T lcm_imp(const T& a, const T& b)
   {
      T temp = boost::math::gcd(a, b);
      return temp ? T(a / temp * b) : T(0);
   }

}

template <typename Integer>
inline BOOST_DEDUCED_TYPENAME enable_if_c<std::numeric_limits<Integer>::is_signed, Integer>::type
lcm(Integer const &a, Integer const &b)
{
   using std::abs;
   return detail::lcm_imp(static_cast<Integer>(abs(a)), static_cast<Integer>(abs(b)));
}


template <typename Integer>
inline BOOST_DEDUCED_TYPENAME disable_if_c<std::numeric_limits<Integer>::is_signed, Integer>::type
lcm(Integer const &a, Integer const &b)
{
   return detail::lcm_imp(a, b);
}


}}
#endif  // BOOST_MATH_COMMON_FACTOR_RT_HPP
