// Copyright John Maddock 2011.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


namespace boost{ namespace math{

template tools::promote_args<BOOST_MATH_TEST_TYPE>::type zeta(BOOST_MATH_TEST_TYPE s, const policies::policy<>&);

template tools::promote_args<BOOST_MATH_TEST_TYPE>::type zeta(BOOST_MATH_TEST_TYPE s);

}} // namespaces
