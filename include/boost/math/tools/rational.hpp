//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_RATIONAL_HPP
#define BOOST_MATH_TOOLS_RATIONAL_HPP

namespace boost{ namespace math{ namespace tools{

namespace detail{

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V&, const mpl::int_<1>*)
{
   return a[0];
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<2>*)
{
   return a[0] + x * a[1];
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<3>*)
{
   return ((a[2] * x) + a[1]) * x + a[0];
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<4>*)
{
   return ((a[3] * x + a[2]) * x + a[1]) * x + a[0];
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<5>*)
{
   return (((a[4] * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0];
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<6>*)
{
   return ((((a[5] * x + a[4]) * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0];
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<7>*)
{
   return (((((a[6] * x + a[5]) * x + a[4]) * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0];
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<8>*)
{
   return ((((((a[7] * x + a[6]) * x + a[5]) * x + a[4]) * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0];
}

template <class T, class V>
inline V evaluate_polynomial_c_imp(const T* a, const V& x, const mpl::int_<9>*)
{
   return (((((((a[8] * x + a[7]) * x + a[6]) * x + a[5]) * x + a[4]) * x + a[3]) * x + a[2]) * x + a[1]) * x + a[0];
}

template <class T, class V, class Tag>
inline V evaluate_polynomial_c_imp(const T* a, const V& val, const Tag*)
{
   return evaluate_polynomial(a, val, Tag::value);
}

} // namespace detail

template <std::size_t N, class T, class V>
inline V evaluate_polynomial(const T(&a)[N], const V& val)
{
   typedef mpl::int_<N> tag_type;
   return detail::evaluate_polynomial_c_imp(static_cast<const T*>(a), val, static_cast<tag_type const*>(0));
}

template <class T, class U>
U evaluate_polynomial(const T* poly, U const& z, std::size_t count)
{
   U sum = static_cast<U>(poly[count - 1]);
   for(int i = static_cast<int>(count) - 2; i >= 0; --i)
   {
      sum *= z;
      sum += static_cast<U>(poly[i]);
   }
   return sum;
}

template <class T, class U>
inline U evaluate_even_polynomial(const T* poly, U z, std::size_t count)
{
   return evaluate_polynomial(poly, z*z, count);
}

template <class T, class U>
inline U evaluate_odd_polynomial(const T* poly, U z, std::size_t count)
{
   return poly[0] + z * evaluate_polynomial(poly+1, z*z, count-1);
}

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

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_RATIONAL_HPP



