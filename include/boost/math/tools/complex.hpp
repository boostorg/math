//  Copyright John Maddock 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//
// Tools for operator on complex as well as scalar types.
//

#ifndef BOOST_MATH_TOOLS_COMPLEX_HPP
#define BOOST_MATH_TOOLS_COMPLEX_HPP

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/is_detected.hpp>

#ifdef BOOST_MATH_ENABLE_CUDA

#include <cuda/std/utility>
#include <cuda/std/complex>

namespace boost {
namespace math {

template <typename T>
using complex = cuda::std::complex<T>;

using cuda::std::real;
using cuda::std::imag;
using cuda::std::abs;
using cuda::std::arg;
using cuda::std::norm;
using cuda::std::conj;
using cuda::std::proj;
using cuda::std::polar;

using cuda::std::exp;
using cuda::std::log;
using cuda::std::log10;

using cuda::std::pow;
using cuda::std::sqrt;

using cuda::std::sin;
using cuda::std::cos;
using cuda::std::tan;
using cuda::std::asin;
using cuda::std::acos;
using cuda::std::atan;

using cuda::std::sinh;
using cuda::std::cosh;
using cuda::std::tanh;
using cuda::std::asinh;
using cuda::std::acosh;
using cuda::std::atanh;

} // namespace math
} // namespace boost

#else

#include <utility>
#include <complex>

namespace boost {
namespace math {

template <typename T>
using complex = std::complex<T>;

using std::real;
using std::imag;
using std::abs;
using std::arg;
using std::norm;
using std::conj;
using std::proj;
using std::polar;

using std::exp;
using std::log;
using std::log10;

using std::pow;
using std::sqrt;

using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;

using std::sinh;
using std::cosh;
using std::tanh;
using std::asinh;
using std::acosh;
using std::atanh;

} // namespace math
} // namespace boost

#endif

namespace boost {
   namespace math {
      namespace tools {

         namespace detail {
         template <typename T, typename = void>
         struct is_complex_type_impl
         {
            static constexpr bool value = false;
         };

         #ifndef BOOST_MATH_ENABLE_CUDA
         template <typename T>
         struct is_complex_type_impl<T, void_t<decltype(std::declval<T>().real()), 
                                               decltype(std::declval<T>().imag())>>
         {
            static constexpr bool value = true;
         };
         #else
         template <typename T>
         struct is_complex_type_impl<T, void_t<decltype(cuda::std::declval<T>().real()), 
                                               decltype(cuda::std::declval<T>().imag())>>
         {
            static constexpr bool value = true;
         };
         #endif
         } // Namespace detail

         template <typename T>
         struct is_complex_type : public detail::is_complex_type_impl<T> {};
         
         //
         // Use this trait to typecast integer literals to something
         // that will interoperate with T:
         //
         template <class T, bool = is_complex_type<T>::value>
         struct integer_scalar_type
         {
            typedef int type;
         };
         template <class T>
         struct integer_scalar_type<T, true>
         {
            typedef typename T::value_type type;
         };
         template <class T, bool = is_complex_type<T>::value>
         struct unsigned_scalar_type
         {
            typedef unsigned type;
         };
         template <class T>
         struct unsigned_scalar_type<T, true>
         {
            typedef typename T::value_type type;
         };
         template <class T, bool = is_complex_type<T>::value>
         struct scalar_type
         {
            typedef T type;
         };
         template <class T>
         struct scalar_type<T, true>
         {
            typedef typename T::value_type type;
         };


} } }

#endif // BOOST_MATH_TOOLS_COMPLEX_HPP
