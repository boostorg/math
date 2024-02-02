//  (C) Copyright Ryan Elandt 2023.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/results_collector.hpp>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/math/concepts/real_concept.hpp>

#ifdef BOOST_HAS_FLOAT128
#include <boost/multiprecision/float128.hpp>
#endif

#if __has_include(<stdfloat>)
#  include <stdfloat>
#endif

using boost::math::tools::detail::ieee754_linear::Layout_IEEE754Linear;
using boost::math::tools::detail::ieee754_linear::Layout_Unspecified;


template <typename T>
bool is_unspeci(Layout_Unspecified<T>) { return true; }

template <typename T>
bool is_paired_with_uint32_t(Layout_IEEE754Linear<T, uint32_t>, T) { return true; }

template <typename T>
bool is_paired_with_uint64_t(Layout_IEEE754Linear<T, uint64_t>, T) { return true; }

#if defined(BOOST_HAS_FLOAT128) && defined(BOOST_HAS_INT128)
template <typename T>
bool is_paired_with_uint128_t(Layout_IEEE754Linear<T, boost::uint128_type>, T) { return true; }
#endif

void test_layout() {
   using boost::math::tools::detail::ieee754_linear::LayoutIdentifier;

   // Check layout for float types
   BOOST_CHECK(is_paired_with_uint32_t(LayoutIdentifier::get_layout<float>(), float()));
   BOOST_CHECK(is_paired_with_uint64_t(LayoutIdentifier::get_layout<double>(), double()));
   BOOST_CHECK(is_unspeci(LayoutIdentifier::get_layout<boost::math::concepts::real_concept>()));
   BOOST_CHECK(is_unspeci(LayoutIdentifier::get_layout<boost::multiprecision::cpp_bin_float_50>()));
   BOOST_CHECK(is_unspeci(LayoutIdentifier::get_layout<boost::multiprecision::cpp_bin_float_100>()));
   BOOST_CHECK(is_unspeci(LayoutIdentifier::get_layout<boost::multiprecision::cpp_dec_float_50>()));
   BOOST_CHECK(is_unspeci(LayoutIdentifier::get_layout<boost::multiprecision::cpp_dec_float_100>()));

// Float128
#if defined(BOOST_HAS_FLOAT128) && defined(BOOST_HAS_INT128)
   BOOST_CHECK(is_paired_with_uint128_t(LayoutIdentifier::get_layout<__float128>(), __float128()));
   BOOST_CHECK(is_paired_with_uint128_t(LayoutIdentifier::get_layout<boost::multiprecision::float128>(), __float128()));
#endif

// stdfloat
#if __has_include(<stdfloat>)
#ifdef __STDCPP_FLOAT32_T__
   BOOST_CHECK(is_paired_with_uint32_t(LayoutIdentifier::get_layout<std::float32_t>(), std::float32_t()));
#endif
#ifdef __STDCPP_FLOAT64_T__
   BOOST_CHECK(is_paired_with_uint64_t(LayoutIdentifier::get_layout<std::float64_t>(), std::float64_t()));
#endif 
#endif

// 80 bit long double
#if LDBL_MANT_DIG == 64 && LDBL_MAX_EXP == 16384
   BOOST_CHECK(is_unspeci(LayoutIdentifier::get_layout<long double>()));
#endif
}

BOOST_AUTO_TEST_CASE( test_main )
{
   test_layout();
   
}
