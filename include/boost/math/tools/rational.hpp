//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_RATIONAL_HPP
#define BOOST_MATH_TOOLS_RATIONAL_HPP

#include <boost/array.hpp>
#include <boost/mpl/int.hpp>

namespace boost{ namespace math{ namespace tools{

//
// Forward declaration to keep two phase lookup happy:
//
template <class T, class U>
U evaluate_polynomial(const T* poly, U const& z, std::size_t count);

namespace detail{
//
// These inline functions evaluate polynomials whose size is
// known at compile time - there's no need for a for-loop
// just an inline application of Horners rule.
//
template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V&, const mpl::int_<0>*)
{
   return static_cast<V>(0);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V&, const mpl::int_<1>*)
{
   return static_cast<V>(a[0]);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<2>*)
{
   return static_cast<V>(a[0] + x * a[1]);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<3>*)
{
   return static_cast<V>(((a[2] * x) + a[1]) * x + a[0]);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<4>*)
{
   return static_cast<V>(((a[3] * x + a[2]) * x + a[1]) * x + a[0]);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<5>*)
{
   return static_cast<V>((((a[4] * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0]);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<6>*)
{
   return static_cast<V>(((((a[5] * x + a[4]) * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0]);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<7>*)
{
   return static_cast<V>((((((a[6] * x + a[5]) * x + a[4]) * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0]);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<8>*)
{
   return static_cast<V>(((((((a[7] * x + a[6]) * x + a[5]) * x + a[4]) * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0]);
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<9>*)
{
   return static_cast<V>((((((((a[8] * x + a[7]) * x + a[6]) * x + a[5]) * x + a[4]) * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0]);
}

template <class T, class V, class Tag>
inline V evaluate_polynomial_c_imp(const T* a, const V& val, const Tag*)
{
   return evaluate_polynomial(a, val, Tag::value);
}

} // namespace detail

//
// Polynomial evaluation with runtime size.
// This requires a for-loop which may be more expensive than
// the loop expanded versions above:
//
template <class T, class U>
U evaluate_polynomial(const T* poly, U const& z, std::size_t count)
{
   BOOST_ASSERT(count > 0);
   U sum = static_cast<U>(poly[count - 1]);
   for(int i = static_cast<int>(count) - 2; i >= 0; --i)
   {
      sum *= z;
      sum += static_cast<U>(poly[i]);
   }
   return sum;
}
//
// Compile time sized polynomials, just inline forwarders to the
// implementations above:
//
template <std::size_t N, class T, class V>
inline V evaluate_polynomial(const T(&a)[N], const V& val)
{
   typedef mpl::int_<N> tag_type;
   return detail::evaluate_polynomial_c_imp(static_cast<const T*>(a), val, static_cast<tag_type const*>(0));
}

template <std::size_t N, class T, class V>
inline V evaluate_polynomial(const boost::array<T,N>& a, const V& val)
{
   typedef mpl::int_<N> tag_type;
   return detail::evaluate_polynomial_c_imp(static_cast<const T*>(a.data()), val, static_cast<tag_type const*>(0));
}
//
// Even polynomials are trivial: just square the argument!
//
template <class T, class U>
inline U evaluate_even_polynomial(const T* poly, U z, std::size_t count)
{
   return evaluate_polynomial(poly, z*z, count);
}

template <std::size_t N, class T, class V>
inline V evaluate_even_polynomial(const T(&a)[N], const V& z)
{
   return evaluate_polynomial(a, z*z);
}

template <std::size_t N, class T, class V>
inline V evaluate_even_polynomial(const boost::array<T,N>& a, const V& z)
{
   return evaluate_polynomial(a, z*z);
}
//
// Odd polynomials come next:
//
template <class T, class U>
inline U evaluate_odd_polynomial(const T* poly, U z, std::size_t count)
{
   return poly[0] + z * evaluate_polynomial(poly+1, z*z, count-1);
}

template <std::size_t N, class T, class V>
inline V evaluate_odd_polynomial(const T(&a)[N], const V& z)
{
   typedef mpl::int_<N-1> tag_type;
   return a[0] + z * detail::evaluate_polynomial_c_imp(static_cast<const T*>(a) + 1, z*z, static_cast<tag_type const*>(0));
}

template <std::size_t N, class T, class V>
inline V evaluate_odd_polynomial(const boost::array<T,N>& a, const V& z)
{
   typedef mpl::int_<N-1> tag_type;
   return a[0] + z * detail::evaluate_polynomial_c_imp(static_cast<const T*>(a.data()) + 1, z*z, static_cast<tag_type const*>(0));
}
//
// Rational functions: numerator and denominator must be
// equal in size.  These always have a for-loop and so may be less
// efficient than evaluating a pair of polynomials. However, there
// are some tricks we can use to prevent overflow that might otherwise
// occur in polynomial evaluation, if z is large.  This is important
// in our Lanczos code for example.
//
template <class T, class U, class V>
V evaluate_rational(const T* num, const U* denom, const V& z_, std::size_t count)
{
   V z(z_);
   V s1, s2;
   if(z <= 1)
   {
      s1 = num[count-1];
      s2 = denom[count-1];
      for(int i = (int)count - 2; i >= 0; --i)
      {
         s1 *= z;
         s2 *= z;
         s1 += num[i];
         s2 += denom[i];
      }
   }
   else
   {
      z = 1 / z;
      s1 = num[0];
      s2 = denom[0];
      for(unsigned i = 1; i < count; ++i)
      {
         s1 *= z;
         s2 *= z;
         s1 += num[i];
         s2 += denom[i];
      }
   }
   return s1 / s2;
}

template <std::size_t N, class T, class U, class V>
inline V evaluate_rational(const T(&a)[N], const U(&b)[N], const V& z)
{
   return evaluate_rational(a, b, z, N);
}

template <std::size_t N, class T, class U, class V>
inline V evaluate_rational(const boost::array<T,N>& a, const boost::array<U,N>& b, const V& z)
{
   return evaluate_rational(a.data(), b.data(), z, N);
}

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_RATIONAL_HPP



