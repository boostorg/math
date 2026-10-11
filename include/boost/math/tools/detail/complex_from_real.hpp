//  Copyright Nick Thompson 2026.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at
//  https://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_DETAIL_COMPLEX_FROM_REAL_HPP
#define BOOST_MATH_TOOLS_DETAIL_COMPLEX_FROM_REAL_HPP

#include <boost/math/tools/config.hpp>

#ifndef BOOST_MATH_BUILD_MODULE
#include <complex>
#include <type_traits>
#endif

// Specialised by each Boost.Multiprecision backend; declared here so that this header does not depend on Multiprecision.
namespace boost { namespace multiprecision {
template <class T>
struct complex_result_from_scalar;
}}

BOOST_MATH_NAMESPACE_BEGIN namespace tools { namespace detail {

// The complex type whose real type is Real, for use with tools::fft and tools::real_fft:
// std::complex<Real> for float, double and long double, and the matching Boost.Multiprecision
// complex type otherwise (cpp_complex_50 for cpp_bin_float_50, and so on). For the latter the
// Multiprecision complex header must be included, e.g. <boost/multiprecision/cpp_complex.hpp>,
// <boost/multiprecision/complex128.hpp> or <boost/multiprecision/mpc.hpp>.
template <class Real, bool = std::is_floating_point<Real>::value>
struct complex_from_real
{
   using type = std::complex<Real>;
};

template <class Real>
struct complex_from_real<Real, false>
{
   using type = typename ::boost::multiprecision::complex_result_from_scalar<Real>::type;
};

}} BOOST_MATH_NAMESPACE_END

#endif
