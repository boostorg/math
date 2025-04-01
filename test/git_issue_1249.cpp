//  (C) Copyright Kilian Kilger 2025.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/test/results_collector.hpp>
#include <boost/math/special_functions/gamma.hpp>

using namespace std;
using namespace boost::math;
using namespace boost::math::policies;

typedef policy<
   policies::domain_error<errno_on_error>,
   policies::pole_error<errno_on_error>,
   policies::overflow_error<errno_on_error>,
   policies::evaluation_error<errno_on_error>
> c_policy;

template<typename T>
struct test_lower
{
   T operator()(T a, T x) const
   {
      return tgamma_lower(a, x, c_policy());
   }

   T expected(T a) const
   {
      return T(0.0);
   }
};

template<typename T>
struct test_upper
{
   T operator()(T a, T x) const
   {
      return tgamma(a, x, c_policy());
   }
   T expected(T a) const
   {
      return tgamma(a, c_policy());
   }
};

template<typename T>
struct test_gamma_p
{
   T operator()(T a, T x) const
   {
      return gamma_p(a, x, c_policy());
   }
   T expected(T) const
   {
      return T(0.0);
   }
};

template<typename T>
struct test_gamma_q
{
   T operator()(T a, T x) const
   {
      return gamma_q(a, x, c_policy());
   }
   T expected(T) const
   {
      return T(1.0);
   }
};

template<typename T, template<typename> class Fun>
void test_impl(T a)
{
   Fun<T> fn;
   errno = 0;
   T x = T(0.0);
   T result = fn(a, x);
   int saveErrno = errno;

   errno = 0;

   T expected = fn.expected(a);

   BOOST_CHECK(errno == saveErrno);
   BOOST_CHECK_EQUAL(result, expected);
}

template<template<typename> class Fun>
void test_type_dispatch(float a)
{
   test_impl<float, Fun>(a);
   test_impl<double, Fun>(double(a));
   test_impl<long double, Fun>(static_cast<long double>(a));
}

template<template<typename> class Fun>
void test_impl()
{
   test_type_dispatch<Fun>(1.0);
   test_type_dispatch<Fun>(0.1);
   test_type_dispatch<Fun>(0.5);
   test_type_dispatch<Fun>(0.6);
   test_type_dispatch<Fun>(1.3);
   test_type_dispatch<Fun>(1.5);
   test_type_dispatch<Fun>(2);
   test_type_dispatch<Fun>(100);
   test_type_dispatch<Fun>(std::numeric_limits<float>::max());
}

BOOST_AUTO_TEST_CASE( test_main )
{
   test_impl<test_lower>();
   test_impl<test_upper>();
   test_impl<test_gamma_p>();
   test_impl<test_gamma_q>();
}