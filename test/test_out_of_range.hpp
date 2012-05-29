// Copyright John Maddock 2012.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TEST_OUT_OF_RANGE_HPP
#define BOOST_MATH_TEST_OUT_OF_RANGE_HPP

#include <boost/math/special_functions/next.hpp>
#include <boost/test/test_tools.hpp>

template <class Distro>
void check_support(const Distro& d)
{
   typedef typename Distro::value_type value_type;
   if((boost::math::isfinite)(range(d).first) && (range(d).first != -boost::math::tools::max_value<value_type>()))
   {
      value_type m = boost::math::float_prior(range(d).first);
      BOOST_ASSERT(m != range(d).first);
      BOOST_ASSERT(m < range(d).first);
      BOOST_CHECK_THROW(pdf(d, m), std::domain_error);
      BOOST_CHECK_THROW(cdf(d, m), std::domain_error);
      BOOST_CHECK_THROW(cdf(complement(d, m)), std::domain_error);
   }
   if((boost::math::isfinite)(range(d).second) && (range(d).second != boost::math::tools::max_value<value_type>()))
   {
      value_type m = boost::math::float_next(range(d).second);
      BOOST_ASSERT(m != range(d).first);
      BOOST_ASSERT(m > range(d).first);
      BOOST_CHECK_THROW(pdf(d, m), std::domain_error);
      BOOST_CHECK_THROW(cdf(d, m), std::domain_error);
      BOOST_CHECK_THROW(cdf(complement(d, m)), std::domain_error);
   }
   if(std::numeric_limits<value_type>::has_infinity)
   {
      if((boost::math::isfinite)(range(d).second))
      {
         BOOST_CHECK_THROW(pdf(d, std::numeric_limits<value_type>::infinity()), std::domain_error);
         BOOST_CHECK_THROW(cdf(d, std::numeric_limits<value_type>::infinity()), std::domain_error);
         BOOST_CHECK_THROW(cdf(complement(d, std::numeric_limits<value_type>::infinity())), std::domain_error);
      }
      if((boost::math::isfinite)(range(d).first))
      {
         BOOST_CHECK_THROW(pdf(d, -std::numeric_limits<value_type>::infinity()), std::domain_error);
         BOOST_CHECK_THROW(cdf(d, -std::numeric_limits<value_type>::infinity()), std::domain_error);
         BOOST_CHECK_THROW(cdf(complement(d, -std::numeric_limits<value_type>::infinity())), std::domain_error);
      }
      BOOST_CHECK_THROW(quantile(d, std::numeric_limits<value_type>::infinity()), std::domain_error);
      BOOST_CHECK_THROW(quantile(d, -std::numeric_limits<value_type>::infinity()), std::domain_error);
      BOOST_CHECK_THROW(quantile(complement(d, std::numeric_limits<value_type>::infinity())), std::domain_error);
      BOOST_CHECK_THROW(quantile(complement(d, -std::numeric_limits<value_type>::infinity())), std::domain_error);
   }
   if(std::numeric_limits<value_type>::has_quiet_NaN)
   {
      BOOST_CHECK_THROW(pdf(d, std::numeric_limits<value_type>::quiet_NaN()), std::domain_error);
      BOOST_CHECK_THROW(cdf(d, std::numeric_limits<value_type>::quiet_NaN()), std::domain_error);
      BOOST_CHECK_THROW(cdf(complement(d, std::numeric_limits<value_type>::quiet_NaN())), std::domain_error);
      BOOST_CHECK_THROW(pdf(d, -std::numeric_limits<value_type>::quiet_NaN()), std::domain_error);
      BOOST_CHECK_THROW(cdf(d, -std::numeric_limits<value_type>::quiet_NaN()), std::domain_error);
      BOOST_CHECK_THROW(cdf(complement(d, -std::numeric_limits<value_type>::quiet_NaN())), std::domain_error);
      BOOST_CHECK_THROW(quantile(d, std::numeric_limits<value_type>::quiet_NaN()), std::domain_error);
      BOOST_CHECK_THROW(quantile(d, -std::numeric_limits<value_type>::quiet_NaN()), std::domain_error);
      BOOST_CHECK_THROW(quantile(complement(d, std::numeric_limits<value_type>::quiet_NaN())), std::domain_error);
      BOOST_CHECK_THROW(quantile(complement(d, -std::numeric_limits<value_type>::quiet_NaN())), std::domain_error);
   }
   BOOST_CHECK_THROW(quantile(d, -1), std::domain_error);
   BOOST_CHECK_THROW(quantile(d, 2), std::domain_error);
   BOOST_CHECK_THROW(quantile(complement(d, -1)), std::domain_error);
   BOOST_CHECK_THROW(quantile(complement(d, 2)), std::domain_error);
}

template <class Distro>
void check_out_of_range()
{
   Distro d;
   check_support(d);
}
template <class Distro>
void check_out_of_range(typename Distro::value_type p1)
{
   typedef typename Distro::value_type value_type;
   Distro d(p1);
   check_support(d);
   if(std::numeric_limits<value_type>::has_infinity)
   {
      BOOST_CHECK_THROW(pdf(Distro(std::numeric_limits<value_type>::infinity()), range(d).first), std::domain_error);
   }
   if(std::numeric_limits<value_type>::has_quiet_NaN)
   {
      BOOST_CHECK_THROW(pdf(Distro(std::numeric_limits<value_type>::quiet_NaN()), range(d).first), std::domain_error);
   }
}
template <class Distro>
void check_out_of_range(typename Distro::value_type p1, typename Distro::value_type p2)
{
   typedef typename Distro::value_type value_type;
   Distro d(p1, p2);
   check_support(d);
   if(std::numeric_limits<value_type>::has_infinity)
   {
      BOOST_CHECK_THROW(pdf(Distro(std::numeric_limits<value_type>::infinity(), p2), range(d).first), std::domain_error);
      BOOST_CHECK_THROW(pdf(Distro(p1, std::numeric_limits<value_type>::infinity()), range(d).first), std::domain_error);
   }
   if(std::numeric_limits<value_type>::has_quiet_NaN)
   {
      BOOST_CHECK_THROW(pdf(Distro(std::numeric_limits<value_type>::quiet_NaN(), p2), range(d).first), std::domain_error);
      BOOST_CHECK_THROW(pdf(Distro(p1, std::numeric_limits<value_type>::quiet_NaN()), range(d).first), std::domain_error);
   }
}

template <class Distro>
void check_out_of_range(typename Distro::value_type p1, typename Distro::value_type p2, typename Distro::value_type p3)
{
   typedef typename Distro::value_type value_type;
   Distro d(p1, p2, p3);
   check_support(d);
   if(std::numeric_limits<value_type>::has_infinity)
   {
      BOOST_CHECK_THROW(pdf(Distro(std::numeric_limits<value_type>::infinity(), p2, p3), range(d).first), std::domain_error);
      BOOST_CHECK_THROW(pdf(Distro(p1, std::numeric_limits<value_type>::infinity(), p3), range(d).first), std::domain_error);
      BOOST_CHECK_THROW(pdf(Distro(p1, p2, std::numeric_limits<value_type>::infinity()), range(d).first), std::domain_error);
   }
   if(std::numeric_limits<value_type>::has_quiet_NaN)
   {
      BOOST_CHECK_THROW(pdf(Distro(std::numeric_limits<value_type>::quiet_NaN(), p2, p3), range(d).first), std::domain_error);
      BOOST_CHECK_THROW(pdf(Distro(p1, std::numeric_limits<value_type>::quiet_NaN(), p3), range(d).first), std::domain_error);
      BOOST_CHECK_THROW(pdf(Distro(p1, p2, std::numeric_limits<value_type>::quiet_NaN()), range(d).first), std::domain_error);
   }
}

#endif // BOOST_MATH_TEST_OUT_OF_RANGE_HPP
