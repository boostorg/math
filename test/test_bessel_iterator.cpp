//  (C) Copyright John Maddock, 2024
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/math/special_functions/bessel_iterators.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>


BOOST_AUTO_TEST_CASE(test_main)
{
   using boost::math::bessel_j_backwards_iterator;
   using boost::math::bessel_i_backwards_iterator;
   using boost::math::bessel_i_forwards_iterator;

   BOOST_CHECK_THROW(bessel_j_backwards_iterator<double>(-2.3, 5.0), std::domain_error);
   BOOST_CHECK_THROW(bessel_j_backwards_iterator<double>(-2.3, 5.0, 2.0), std::domain_error);
   BOOST_CHECK_THROW(bessel_j_backwards_iterator<double>(-2.3, 5.0, 2.0, 2.0), std::domain_error);

   double tolerance = std::numeric_limits<double>::epsilon() * 5;

   {
      bessel_j_backwards_iterator<double> i(12, 2.0);
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_j(8, 2.0), tolerance);
      double v = *i++;
      BOOST_CHECK_CLOSE_FRACTION(v, boost::math::cyl_bessel_j(8, 2.0), tolerance);
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_j(7, 2.0), tolerance);
   }
   {
      bessel_j_backwards_iterator<double> i(12, 2.0, boost::math::cyl_bessel_j(12, 2.0));
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_j(8, 2.0), tolerance);
   }
   {
      bessel_j_backwards_iterator<double> i(12, 2.0, boost::math::cyl_bessel_j(13, 2.0), boost::math::cyl_bessel_j(12, 2.0));
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_j(8, 2.0), tolerance);
   }

   {
      bessel_i_backwards_iterator<double> i(12, 2.0);
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(8, 2.0), tolerance);
      double v = *i++;
      BOOST_CHECK_CLOSE_FRACTION(v, boost::math::cyl_bessel_i(8, 2.0), tolerance);
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(7, 2.0), tolerance);

      BOOST_CHECK_THROW(bessel_i_backwards_iterator<double>(-1.1, 2.0), std::domain_error);
      BOOST_CHECK_THROW(bessel_i_backwards_iterator<double>(-1.1, 2.0, 1.0), std::domain_error);
      BOOST_CHECK_THROW(bessel_i_backwards_iterator<double>(-1.1, 2.0, 1.0, 1.0), std::domain_error);
   }
   {
      bessel_i_backwards_iterator<double> i(12, 2.0, boost::math::cyl_bessel_i(12, 2.0));
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(8, 2.0), tolerance);
   }
   {
      bessel_i_backwards_iterator<double> i(12, 2.0, boost::math::cyl_bessel_i(13, 2.0), boost::math::cyl_bessel_i(12, 2.0));
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(8, 2.0), tolerance);
   }

   {
      bessel_i_forwards_iterator<double> i(-12, 2.0);
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(-8, 2.0), tolerance);
      double v = *i++;
      BOOST_CHECK_CLOSE_FRACTION(v, boost::math::cyl_bessel_i(-8, 2.0), tolerance);
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(-7, 2.0), tolerance);

      BOOST_CHECK_THROW(bessel_i_forwards_iterator<double>(1.1, 2.0), std::domain_error);
      BOOST_CHECK_THROW(bessel_i_forwards_iterator<double>(1.1, 2.0, 1.0), std::domain_error);
      BOOST_CHECK_THROW(bessel_i_forwards_iterator<double>(1.1, 2.0, 1.0, 1.0), std::domain_error);
   }
   {
      bessel_i_forwards_iterator<double> i(-12, 2.0, boost::math::cyl_bessel_i(-12, 2.0));
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(-8, 2.0), tolerance);
   }
   {
      bessel_i_forwards_iterator<double> i(-12, 2.0, boost::math::cyl_bessel_i(-13, 2.0), boost::math::cyl_bessel_i(-12, 2.0));
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(-8, 2.0), tolerance);
   }
}

BOOST_AUTO_TEST_CASE(test_order_offset)
{
   using boost::math::bessel_j_backwards_iterator;
   using boost::math::bessel_i_backwards_iterator;
   using boost::math::bessel_i_forwards_iterator;
   using boost::math::tools::recurrence_offset;
   using mp_t = boost::multiprecision::cpp_bin_float_50;
   // I_{v+j}(x) and J_{v+j}(x) to double precision, with v + j exact:
   auto I = [](double v, int j, double x) { return static_cast<double>(boost::math::cyl_bessel_i(mp_t(v) + j, mp_t(x))); };
   auto J = [](double v, int j, double x) { return static_cast<double>(boost::math::cyl_bessel_j(mp_t(v) + j, mp_t(x))); };

   //
   // Starting from order v + k with k an integer, each order v + j is rounded only once.  Starting from the rounded
   // order v + k instead puts every order out by ~eps * k; for I near v = -1 and small x that costs over 5000 eps:
   //
   double tolerance = std::numeric_limits<double>::epsilon() * 20;
   {
      double v = -0.999, x = 0.01;
      bessel_i_backwards_iterator<double> i(v, x, I(v, 21, x), I(v, 20, x), recurrence_offset(20));
      for (int j = 20; j >= 0; --j, ++i)
         BOOST_CHECK_CLOSE_FRACTION(*i, I(v, j, x), tolerance);
      bessel_i_backwards_iterator<double> i1(v, x, I(v, 20, x), recurrence_offset(20));
      for (int j = 20; j >= 0; --j, ++i1)
         BOOST_CHECK_CLOSE_FRACTION(*i1, I(v, j, x), tolerance);
   }
   {
      double v = 0.1, x = 2;
      bessel_i_backwards_iterator<double> i(v, x, I(v, 101, x), I(v, 100, x), recurrence_offset(100));
      for (int j = 100; j >= 0; --j, ++i)
         BOOST_CHECK_CLOSE_FRACTION(*i, I(v, j, x), tolerance);
   }
   {
      double v = 0.1, x = 2;
      bessel_j_backwards_iterator<double> i(v, x, J(v, 41, x), J(v, 40, x), recurrence_offset(40));
      for (int j = 40; j >= 0; --j, ++i)
         BOOST_CHECK_CLOSE_FRACTION(*i, J(v, j, x), tolerance);
      bessel_j_backwards_iterator<double> i1(v, x, J(v, 40, x), recurrence_offset(40));
      for (int j = 40; j >= 0; --j, ++i1)
         BOOST_CHECK_CLOSE_FRACTION(*i1, J(v, j, x), tolerance);
   }
   //
   // The constructors that compute the starting value, and the forwards iterator (only stable for integer orders):
   //
   tolerance = std::numeric_limits<double>::epsilon() * 5;
   {
      bessel_i_backwards_iterator<double> i(2.0, 2.0, recurrence_offset(10));
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(8, 2.0), tolerance);
      bessel_j_backwards_iterator<double> k(2.0, 2.0, recurrence_offset(10));
      for (unsigned j = 0; j < 4; ++j)
         ++k;
      BOOST_CHECK_CLOSE_FRACTION(*k, boost::math::cyl_bessel_j(8, 2.0), tolerance);
   }
   {
      bessel_i_forwards_iterator<double> i(0.0, 2.0, recurrence_offset(-12));
      for (unsigned j = 0; j < 4; ++j)
         ++i;
      BOOST_CHECK_CLOSE_FRACTION(*i, boost::math::cyl_bessel_i(-8, 2.0), tolerance);
      bessel_i_forwards_iterator<double> i1(0.0, 2.0, boost::math::cyl_bessel_i(-12, 2.0), recurrence_offset(-12));
      for (unsigned j = 0; j < 4; ++j)
         ++i1;
      BOOST_CHECK_CLOSE_FRACTION(*i1, boost::math::cyl_bessel_i(-8, 2.0), tolerance);
      bessel_i_forwards_iterator<double> i2(0.0, 2.0, boost::math::cyl_bessel_i(-13, 2.0), boost::math::cyl_bessel_i(-12, 2.0), recurrence_offset(-12));
      for (unsigned j = 0; j < 4; ++j)
         ++i2;
      BOOST_CHECK_CLOSE_FRACTION(*i2, boost::math::cyl_bessel_i(-8, 2.0), tolerance);
   }
   // The domain checks apply to the starting order v + k:
   BOOST_CHECK_THROW(bessel_i_backwards_iterator<double>(0.5, 2.0, 1.0, recurrence_offset(-2)), std::domain_error);
   BOOST_CHECK_THROW(bessel_i_backwards_iterator<double>(0.5, 2.0, 1.0, 1.0, recurrence_offset(-2)), std::domain_error);
   BOOST_CHECK_THROW(bessel_j_backwards_iterator<double>(0.5, 2.0, 1.0, recurrence_offset(-1)), std::domain_error);
   BOOST_CHECK_THROW(bessel_i_forwards_iterator<double>(-0.5, 2.0, 1.0, recurrence_offset(2)), std::domain_error);
}
