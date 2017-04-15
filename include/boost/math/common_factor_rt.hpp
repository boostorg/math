//  (C) Copyright Jeremy William Murphy 2016.

//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_COMMON_FACTOR_RT_HPP
#define BOOST_MATH_COMMON_FACTOR_RT_HPP

#include <boost/assert.hpp>
#include <boost/core/enable_if.hpp>

#include <boost/config.hpp>  // for BOOST_NESTED_TEMPLATE, etc.
#include <boost/limits.hpp>  // for std::numeric_limits
#include <climits>           // for CHAR_MIN
#include <boost/detail/workaround.hpp>
#include <iterator>
#include <algorithm>
#include <limits>
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
#include <type_traits>
#endif

#if (defined(BOOST_MSVC) || (defined(__clang__) && defined(__c2__)) || (defined(BOOST_INTEL) && defined(_MSC_VER))) && (defined(_M_IX86) || defined(_M_X64))
#include <intrin.h>
#endif

#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable:4127 4244)  // Conditional expression is constant
#endif

#if !defined(BOOST_NO_CXX11_HDR_TYPE_TRAITS) && !defined(BOOST_NO_CXX11_NOEXCEPT)
#define BOOST_GCD_NOEXCEPT(T) noexcept(std::is_arithmetic<T>::value)
#else
#define BOOST_GCD_NOEXCEPT(T)
#endif

namespace boost {

   template <class I>
   class rational;

   namespace math {

      namespace gcd_detail{

         //
         // some helper functions which really should be constexpr already, but sadly aren't:
         //
#ifndef BOOST_NO_CXX14_CONSTEXPR
         template <class T>
         inline constexpr auto constexpr_swap(T&a, T& b) BOOST_GCD_NOEXCEPT(T) -> decltype(a.swap(b))
         {
            return a.swap(b);
         }
         template <class T, class U>
         inline constexpr void constexpr_swap(T&a, U& b...) BOOST_GCD_NOEXCEPT(T)
         {
            T t(static_cast<T&&>(a));
            a = static_cast<T&&>(b);
            b = static_cast<T&&>(t);
         }
#else
         template <class T>
         inline void constexpr_swap(T&a, T& b) BOOST_GCD_NOEXCEPT(T)
         {
            using std::swap;
            swap(a, b);
         }
#endif

      template <class T, bool a = 
#ifndef BOOST_NO_CXX11_HDR_TYPE_TRAITS
         std::is_unsigned<T>::value || 
#endif
         (std::numeric_limits<T>::is_specialized && !std::numeric_limits<T>::is_signed)>
      struct gcd_traits_abs_defaults
      {
         inline static BOOST_CXX14_CONSTEXPR const T& abs(const T& val) BOOST_GCD_NOEXCEPT(T) { return val; }
      };
      template <class T>
      struct gcd_traits_abs_defaults<T, false>
      {
         inline static T BOOST_CXX14_CONSTEXPR abs(const T& val) BOOST_GCD_NOEXCEPT(T)
         {
            // This sucks, but std::abs is not constexpr :(
            return val < 0 ? -val : val;
         }
      };

      enum method_type
      {
         method_euclid = 0,
         method_binary = 1,
         method_mixed = 2,
      };


      template <class T>
      struct gcd_traits_defaults : public gcd_traits_abs_defaults<T>
      {
         BOOST_FORCEINLINE static BOOST_CXX14_CONSTEXPR unsigned make_odd(T& val) BOOST_GCD_NOEXCEPT(T)
         {
            unsigned r = 0;
            while(0 == (val & 1u))
            {
#ifdef _MSC_VER  // VC++ can't handle operator >>= in constexpr code for some reason
               val = val >> 1;
#else
               val >>= 1;
#endif
               ++r;
            }
            return r;
         }
         inline static BOOST_CXX14_CONSTEXPR bool less(const T& a, const T& b) BOOST_GCD_NOEXCEPT(T)
         {
            return a < b;
         }
         static const method_type method = std::numeric_limits<T>::is_specialized && std::numeric_limits<T>::is_integer ? method_mixed : method_euclid;
      };
      //
      // Default gcd_traits just inherits from defaults:
      //
      template <class T>
      struct gcd_traits : public gcd_traits_defaults<T> {};

#ifdef BOOST_NO_CXX14_CONSTEXPR
      //
      // Some platforms have fast bitscan operations, that allow us to implement
      // make_odd much more efficiently, unfortunately we can't use these if we want
      // the functions to be constexpr as the compiler intrinsics aren't constexpr.
      //
#if (defined(BOOST_MSVC) || (defined(__clang__) && defined(__c2__)) || (defined(BOOST_INTEL) && defined(_MSC_VER))) && (defined(_M_IX86) || defined(_M_X64))
#pragma intrinsic(_BitScanForward,)
      template <>
      struct gcd_traits<unsigned long> : public gcd_traits_defaults<unsigned long>
      {
         BOOST_FORCEINLINE static unsigned find_lsb(unsigned long val) BOOST_NOEXCEPT
         {
            unsigned long result;
            _BitScanForward(&result, val);
            return result;
         }
         BOOST_FORCEINLINE static unsigned make_odd(unsigned long& val) BOOST_NOEXCEPT
         {
            unsigned result = find_lsb(val);
            val >>= result;
            return result;
         }
      };

#ifdef _M_X64
#pragma intrinsic(_BitScanForward64)
      template <>
      struct gcd_traits<unsigned __int64> : public gcd_traits_defaults<unsigned __int64>
      {
         BOOST_FORCEINLINE static unsigned find_lsb(unsigned __int64 mask) BOOST_NOEXCEPT
         {
            unsigned long result;
            _BitScanForward64(&result, mask);
            return result;
         }
         BOOST_FORCEINLINE static unsigned make_odd(unsigned __int64& val) BOOST_NOEXCEPT
         {
            unsigned result = find_lsb(val);
            val >>= result;
            return result;
         }
      };
#endif
      //
      // Other integer type are trivial adaptations of the above,
      // this works for signed types too, as by the time these functions
      // are called, all values are > 0.
      //
      template <> struct gcd_traits<long> : public gcd_traits_defaults<long> 
      { BOOST_FORCEINLINE static unsigned make_odd(long& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
      template <> struct gcd_traits<unsigned int> : public gcd_traits_defaults<unsigned int> 
      { BOOST_FORCEINLINE static unsigned make_odd(unsigned int& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
      template <> struct gcd_traits<int> : public gcd_traits_defaults<int> 
      { BOOST_FORCEINLINE static unsigned make_odd(int& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
      template <> struct gcd_traits<unsigned short> : public gcd_traits_defaults<unsigned short> 
      { BOOST_FORCEINLINE static unsigned make_odd(unsigned short& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
      template <> struct gcd_traits<short> : public gcd_traits_defaults<short> 
      { BOOST_FORCEINLINE static unsigned make_odd(short& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
      template <> struct gcd_traits<unsigned char> : public gcd_traits_defaults<unsigned char> 
      { BOOST_FORCEINLINE static unsigned make_odd(unsigned char& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
      template <> struct gcd_traits<signed char> : public gcd_traits_defaults<signed char> 
      { BOOST_FORCEINLINE static signed make_odd(signed char& val)BOOST_NOEXCEPT{ signed result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
      template <> struct gcd_traits<char> : public gcd_traits_defaults<char> 
      { BOOST_FORCEINLINE static unsigned make_odd(char& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
      template <> struct gcd_traits<wchar_t> : public gcd_traits_defaults<wchar_t> 
      { BOOST_FORCEINLINE static unsigned make_odd(wchar_t& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; } };
#endif
#ifdef _M_X64
      template <> struct gcd_traits<__int64> : public gcd_traits_defaults<__int64> 
      { BOOST_FORCEINLINE static unsigned make_odd(__int64& val)BOOST_NOEXCEPT{ unsigned result = gcd_traits<unsigned __int64>::find_lsb(val); val >>= result; return result; } };
#endif

#elif defined(BOOST_GCC) || defined(__clang__) || (defined(BOOST_INTEL) && defined(__GNUC__))

      template <>
      struct gcd_traits<unsigned> : public gcd_traits_defaults<unsigned>
      {
         BOOST_FORCEINLINE static unsigned find_lsb(unsigned mask)BOOST_NOEXCEPT
         {
            return __builtin_ctz(mask);
         }
         BOOST_FORCEINLINE static unsigned make_odd(unsigned& val)BOOST_NOEXCEPT
         {
            unsigned result = find_lsb(val);
            val >>= result;
            return result;
         }
      };
      template <>
      struct gcd_traits<unsigned long> : public gcd_traits_defaults<unsigned long>
      {
         BOOST_FORCEINLINE static unsigned find_lsb(unsigned long mask)BOOST_NOEXCEPT
         {
            return __builtin_ctzl(mask);
         }
         BOOST_FORCEINLINE static unsigned make_odd(unsigned long& val)BOOST_NOEXCEPT
         {
            unsigned result = find_lsb(val);
            val >>= result;
            return result;
         }
      };
      template <>
      struct gcd_traits<boost::ulong_long_type> : public gcd_traits_defaults<boost::ulong_long_type>
      {
         BOOST_FORCEINLINE static unsigned find_lsb(boost::ulong_long_type mask)BOOST_NOEXCEPT
         {
            return __builtin_ctzll(mask);
         }
         BOOST_FORCEINLINE static unsigned make_odd(boost::ulong_long_type& val)BOOST_NOEXCEPT
         {
            unsigned result = find_lsb(val);
            val >>= result;
            return result;
         }
      };
      //
      // Other integer type are trivial adaptations of the above,
      // this works for signed types too, as by the time these functions
      // are called, all values are > 0.
      //
      template <> struct gcd_traits<boost::long_long_type> : public gcd_traits_defaults<boost::long_long_type>
      {
         BOOST_FORCEINLINE static unsigned make_odd(boost::long_long_type& val)BOOST_NOEXCEPT { unsigned result = gcd_traits<boost::ulong_long_type>::find_lsb(val); val >>= result; return result; }
      };
      template <> struct gcd_traits<long> : public gcd_traits_defaults<long>
      {
         BOOST_FORCEINLINE static unsigned make_odd(long& val)BOOST_NOEXCEPT { unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; }
      };
      template <> struct gcd_traits<int> : public gcd_traits_defaults<int>
      {
         BOOST_FORCEINLINE static unsigned make_odd(int& val)BOOST_NOEXCEPT { unsigned result = gcd_traits<unsigned long>::find_lsb(val); val >>= result; return result; }
      };
      template <> struct gcd_traits<unsigned short> : public gcd_traits_defaults<unsigned short>
      {
         BOOST_FORCEINLINE static unsigned make_odd(unsigned short& val)BOOST_NOEXCEPT { unsigned result = gcd_traits<unsigned>::find_lsb(val); val >>= result; return result; }
      };
      template <> struct gcd_traits<short> : public gcd_traits_defaults<short>
      {
         BOOST_FORCEINLINE static unsigned make_odd(short& val)BOOST_NOEXCEPT { unsigned result = gcd_traits<unsigned>::find_lsb(val); val >>= result; return result; }
      };
      template <> struct gcd_traits<unsigned char> : public gcd_traits_defaults<unsigned char>
      {
         BOOST_FORCEINLINE static unsigned make_odd(unsigned char& val)BOOST_NOEXCEPT { unsigned result = gcd_traits<unsigned>::find_lsb(val); val >>= result; return result; }
      };
      template <> struct gcd_traits<signed char> : public gcd_traits_defaults<signed char>
      {
         BOOST_FORCEINLINE static signed make_odd(signed char& val)BOOST_NOEXCEPT { signed result = gcd_traits<unsigned>::find_lsb(val); val >>= result; return result; }
      };
      template <> struct gcd_traits<char> : public gcd_traits_defaults<char>
      {
         BOOST_FORCEINLINE static unsigned make_odd(char& val)BOOST_NOEXCEPT { unsigned result = gcd_traits<unsigned>::find_lsb(val); val >>= result; return result; }
      };
#ifndef BOOST_NO_INTRINSIC_WCHAR_T
      template <> struct gcd_traits<wchar_t> : public gcd_traits_defaults<wchar_t>
      {
         BOOST_FORCEINLINE static unsigned make_odd(wchar_t& val)BOOST_NOEXCEPT { unsigned result = gcd_traits<unsigned>::find_lsb(val); val >>= result; return result; }
      };
#endif
#endif
#endif // BOOST_NO_CXX14_CONSTEXPR
   //
   // The Mixed Binary Euclid Algorithm
   // Sidi Mohamed Sedjelmaci
   // Electronic Notes in Discrete Mathematics 35 (2009) 169-176
   //
   template <class T>
   BOOST_CXX14_CONSTEXPR T mixed_binary_gcd(T u, T v) BOOST_GCD_NOEXCEPT(T)
   {
      if(gcd_traits<T>::less(u, v))
         constexpr_swap(u, v);

      unsigned shifts = 0;

      if(!u)
         return v;
      if(!v)
         return u;

      shifts = (std::min)(gcd_traits<T>::make_odd(u), gcd_traits<T>::make_odd(v));

      while(gcd_traits<T>::less(1, v))
      {
         u %= v;
         v -= u;
         if(!u)
            return v << shifts;
         if(!v)
            return u << shifts;
         gcd_traits<T>::make_odd(u);
         gcd_traits<T>::make_odd(v);
         if(gcd_traits<T>::less(u, v))
            constexpr_swap(u, v);
      }
      return (v == 1 ? v : u) << shifts;
   }

    /** Stein gcd (aka 'binary gcd')
     * 
     * From Mathematics to Generic Programming, Alexander Stepanov, Daniel Rose
     */
    template <typename SteinDomain>
    BOOST_CXX14_CONSTEXPR SteinDomain Stein_gcd(SteinDomain m, SteinDomain n) BOOST_GCD_NOEXCEPT(SteinDomain)
    {
        BOOST_ASSERT(m >= 0);
        BOOST_ASSERT(n >= 0);
        if (m == SteinDomain(0))
            return n;
        if (n == SteinDomain(0))
            return m;
        // m > 0 && n > 0
        int d_m = gcd_traits<SteinDomain>::make_odd(m);
        int d_n = gcd_traits<SteinDomain>::make_odd(n);
        // odd(m) && odd(n)
        while (m != n)
        {
            if (n > m)
               constexpr_swap(n, m);
            m -= n;
            gcd_traits<SteinDomain>::make_odd(m);
        }
        // m == n
        m <<= (std::min)(d_m, d_n);
        return m;
    }

    
    /** Euclidean algorithm
     * 
     * From Mathematics to Generic Programming, Alexander Stepanov, Daniel Rose
     * 
     */
    template <typename EuclideanDomain>
    inline BOOST_CXX14_CONSTEXPR EuclideanDomain Euclid_gcd(EuclideanDomain a, EuclideanDomain b) BOOST_GCD_NOEXCEPT(EuclideanDomain)
    {
        while (b != EuclideanDomain(0))
        {
            a %= b;
            constexpr_swap(a, b);
        }
        return a;
    }


    template <typename T>
    inline BOOST_CXX14_CONSTEXPR BOOST_DEDUCED_TYPENAME enable_if_c<gcd_traits<T>::method == method_mixed, T>::type
       optimal_gcd_select(T const &a, T const &b) BOOST_GCD_NOEXCEPT(T)
    {
       return gcd_detail::mixed_binary_gcd(a, b);
    }

    template <typename T>
    inline BOOST_CXX14_CONSTEXPR BOOST_DEDUCED_TYPENAME enable_if_c<gcd_traits<T>::method == method_binary, T>::type
       optimal_gcd_select(T const &a, T const &b) BOOST_GCD_NOEXCEPT(T)
    {
       return gcd_detail::Stein_gcd(a, b);
    }

    template <typename T>
    inline BOOST_CXX14_CONSTEXPR BOOST_DEDUCED_TYPENAME enable_if_c<gcd_traits<T>::method == method_euclid, T>::type
       optimal_gcd_select(T const &a, T const &b) BOOST_GCD_NOEXCEPT(T)
    {
       return gcd_detail::Euclid_gcd(a, b);
    }

    template <class T>
    inline BOOST_CXX14_CONSTEXPR T lcm_imp(const T& a, const T& b) BOOST_GCD_NOEXCEPT(T)
    {
       T temp = boost::math::gcd_detail::optimal_gcd_select(a, b);
#if BOOST_WORKAROUND(BOOST_GCC_VERSION, < 40500)
       return (temp != T(0)) ? T(a / temp * b) : T(0);
#else
       return temp ? T(a / temp * b) : T(0);
#endif
    }

} // namespace detail


template <typename Integer>
inline BOOST_CXX14_CONSTEXPR Integer gcd(Integer const &a, Integer const &b) BOOST_GCD_NOEXCEPT(Integer)
{
    return gcd_detail::optimal_gcd_select(static_cast<Integer>(gcd_detail::gcd_traits<Integer>::abs(a)), static_cast<Integer>(gcd_detail::gcd_traits<Integer>::abs(b)));
}

template <typename Integer>
inline BOOST_CXX14_CONSTEXPR Integer lcm(Integer const &a, Integer const &b) BOOST_GCD_NOEXCEPT(Integer)
{
   return gcd_detail::lcm_imp(static_cast<Integer>(gcd_detail::gcd_traits<Integer>::abs(a)), static_cast<Integer>(gcd_detail::gcd_traits<Integer>::abs(b)));
}

template <typename Integer, typename... Args>
inline BOOST_CXX14_CONSTEXPR Integer gcd(Integer const &a, Integer const &b, Args const&... args) BOOST_GCD_NOEXCEPT(Integer)
{
   return gcd(a, gcd(b, args...));
}

template <typename Integer, typename... Args>
inline BOOST_CXX14_CONSTEXPR Integer lcm(Integer const &a, Integer const &b, Args const&... args) BOOST_GCD_NOEXCEPT(Integer)
{
   return lcm(a, lcm(b, args...));
}

//
// Special handling for rationals:
//
template <typename Integer>
inline typename boost::enable_if_c<std::numeric_limits<Integer>::is_specialized, boost::rational<Integer> >::type gcd(boost::rational<Integer> const &a, boost::rational<Integer> const &b)
{
   return boost::rational<Integer>(static_cast<Integer>(gcd(a.numerator(), b.numerator())), static_cast<Integer>(lcm(a.denominator(), b.denominator())));
}

template <typename Integer>
inline typename boost::enable_if_c<std::numeric_limits<Integer>::is_specialized, boost::rational<Integer> >::type lcm(boost::rational<Integer> const &a, boost::rational<Integer> const &b)
{
   return boost::rational<Integer>(static_cast<Integer>(lcm(a.numerator(), b.numerator())), static_cast<Integer>(gcd(a.denominator(), b.denominator())));
}
/**
 * Knuth, The Art of Computer Programming: Volume 2, Third edition, 1998
 * Chapter 4.5.2, Algorithm C: Greatest common divisor of n integers.
 *
 * Knuth counts down from n to zero but we naturally go from first to last.
 * We also return the termination position because it might be useful to know.
 * 
 * Partly by quirk, partly by design, this algorithm is defined for n = 1, 
 * because the gcd of {x} is x. It is not defined for n = 0.
 * 
 * @tparam  I   Input iterator.
 * @return  The gcd of the range and the iterator position at termination.
 */
template <typename I>
std::pair<typename std::iterator_traits<I>::value_type, I>
gcd_range(I first, I last) BOOST_GCD_NOEXCEPT(I)
{
    BOOST_ASSERT(first != last);
    typedef typename std::iterator_traits<I>::value_type T;
    
    T d = *first++;
    while (d != T(1) && first != last)
    {
        d = gcd(d, *first);
        first++;
    }
    return std::make_pair(d, first);
}
template <typename I>
std::pair<typename std::iterator_traits<I>::value_type, I>
lcm_range(I first, I last) BOOST_GCD_NOEXCEPT(I)
{
    BOOST_ASSERT(first != last);
    typedef typename std::iterator_traits<I>::value_type T;
    
    T d = *first++;
    while (d != T(1) && first != last)
    {
        d = lcm(d, *first);
        first++;
    }
    return std::make_pair(d, first);
}

template < typename IntegerType >
struct gcd_evaluator
#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL
   : public std::binary_function<IntegerType, IntegerType, IntegerType>
#endif
{
#ifndef BOOST_NO_CXX11_HDR_FUNCTIONAL
   typedef IntegerType first_argument_type;
   typedef IntegerType second_argument_type;
   typedef IntegerType result_type;
#endif
   IntegerType operator()(IntegerType const &a, IntegerType const &b)const
   {
      return boost::math::gcd(a, b);
   }
};

template < typename IntegerType >
struct lcm_evaluator
#ifdef BOOST_NO_CXX11_HDR_FUNCTIONAL
   : public std::binary_function<IntegerType, IntegerType, IntegerType>
#endif
{
#ifndef BOOST_NO_CXX11_HDR_FUNCTIONAL
   typedef IntegerType first_argument_type;
   typedef IntegerType second_argument_type;
   typedef IntegerType result_type;
#endif
   IntegerType operator()(IntegerType const &a, IntegerType const &b)const
   {
      return boost::math::lcm(a, b);
   }
};

}  // namespace math
}  // namespace boost

#ifdef BOOST_MSVC
#pragma warning(pop)
#endif

#endif  // BOOST_MATH_COMMON_FACTOR_RT_HPP
