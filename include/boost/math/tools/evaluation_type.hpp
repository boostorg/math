//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_EVALUATION_TYPE_HPP
#define BOOST_MATH_TOOLS_EVALUATION_TYPE_HPP

namespace boost{ namespace math{ namespace tools{

template <class T>
struct evaluation
{
   typedef T type;
};

#ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
#define BOOST_MATH_FLT_DIGITS std::numeric_limits<float>::digits
#define BOOST_MATH_DBL_DIGITS std::numeric_limits<double>::digits
#define BOOST_MATH_LDBL_DIGITS std::numeric_limits<long double>::digits
#else
#define BOOST_MATH_FLT_DIGITS FLT_MANT_DIG
#define BOOST_MATH_DBL_DIGITS DBL_MANT_DIG
#define BOOST_MATH_LDBL_DIGITS LDBL_MANT_DIG
#endif

template<>
struct evaluation<float>
{
   typedef boost::mpl::if_c<
      BOOST_MATH_DBL_DIGITS >= 40,
      double, float>::type type;
};

template<>
struct evaluation<double>
{
   typedef boost::mpl::if_c<
      BOOST_MATH_LDBL_DIGITS >= 64,
      long double,
      double>::type type;
};

}}} // namespaces

#endif // BOOST_MATH_TOOLS_EVALUATION_TYPE_HPP



