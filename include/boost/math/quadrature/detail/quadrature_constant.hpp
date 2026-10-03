//  Copyright John Maddock 2015.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


#ifndef BOOST_MATH_QUADRATURE_DETAIL_QUADRATURE_CONSTANT_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_QUADRATURE_CONSTANT_HPP

#ifndef BOOST_MATH_BUILD_MODULE
#include <limits>
#include <type_traits>
#endif

namespace boost { namespace math{ namespace quadrature{ namespace detail {

template <class T>
struct quadrature_constant_category
{
   static const unsigned value =
      (std::numeric_limits<T>::is_specialized == 0) ? 999 :
      (std::numeric_limits<T>::radix == 2) ?
      (
#ifdef BOOST_HAS_FLOAT128
         (std::numeric_limits<T>::digits <= 113) && std::is_constructible<T, __float128>::value ? 0 :
#else
         (std::numeric_limits<T>::digits <= std::numeric_limits<long double>::digits) && std::is_constructible<T, long double>::value ? 0 :
#endif
         (std::numeric_limits<T>::digits10 <= 110) && std::is_constructible<T, const char*>::value ? 4 : 999
      ) : (std::numeric_limits<T>::digits10 <= 110) && std::is_constructible<T, const char*>::value ? 4 : 999;

   using storage_type =
      std::conditional_t<(std::numeric_limits<T>::is_specialized == 0), T,
         std::conditional_t<(std::numeric_limits<T>::radix == 2),
            std::conditional_t< ((std::numeric_limits<T>::digits <= std::numeric_limits<float>::digits) && std::is_constructible<T, float>::value),
               float,
               std::conditional_t<((std::numeric_limits<T>::digits <= std::numeric_limits<double>::digits) && std::is_constructible<T, double>::value),
                  double,
                  std::conditional_t<((std::numeric_limits<T>::digits <= std::numeric_limits<long double>::digits) && std::is_constructible<T, long double>::value),
                     long double,
#ifdef BOOST_HAS_FLOAT128
                     std::conditional_t<((std::numeric_limits<T>::digits <= 113) && std::is_constructible<T, __float128>::value),
                        __float128,
                        T
                     >
                  >
#else
                     T
                  >
#endif
               >
            >, T
         >
      >;
};

}
}
}
}

#endif // BOOST_MATH_QUADRATURE_DETAIL_QUADRATURE_CONSTANT_HPP