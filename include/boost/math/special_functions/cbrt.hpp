//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_CBRT_HPP
#define BOOST_MATH_SF_CBRT_HPP

#include <boost/math/tools/roots.hpp>

namespace boost{ namespace math{

namespace detail
{

	template <class T>
	struct cbrt_functor
	{
		 cbrt_functor(T const& target) : a(target){}
		 std::tr1::tuple<T, T, T> operator()(T const& z)
		 {
				T sqr = z * z;
				return std::tr1::make_tuple(sqr * z - a, 3 * sqr, 6 * z);
		 }
	private:
		 T a;
	};

} // namespace detail

template <class T>
T cbrt(T z)
{
   using namespace std;
   int exp, sign(1);
   if(z < 0)
   {
      z = -z;
      sign = -sign;
   }
   if(z == 0)
      return 0;

   frexp(z, &exp);
   T min = static_cast<T>(ldexp(0.5, exp/3));
   T max = static_cast<T>(ldexp(2.0, exp/3));
   T guess = static_cast<T>(ldexp(1.0, exp/3));
   int digits = (tools::digits<T>()) / 2;
   return sign * tools::halley_iterate(detail::cbrt_functor<T>(z), guess, min, max, digits);
}

} // namespace math
} // namespace boost

#endif // BOOST_MATH_SF_CBRT_HPP


