//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_TEST_HPP
#define BOOST_MATH_TOOLS_TEST_HPP

#include <boost/math/tools/stats.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#if defined(__CYGWIN__)
#define BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
#endif

namespace boost{ namespace math{ namespace tools{

template <class T>
struct test_result
{
   boost::math::tools::stats<T> stat;   // statistics for the test
   unsigned worst_case;                 // index of the worst case test

   test_result& operator+=(const test_result& t)
   {
      if((t.stat.max)() > (stat.max)())
         worst_case = t.worst_case;
      stat += t.stat;
      return *this;
   }
};

template <class T>
struct calculate_result_type
{
   typedef typename T::value_type row_type;
   typedef typename row_type::value_type value_type;
};

template <class T>
T relative_error(T a, T b)
{
   using namespace std;
#ifdef BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS
   //
   // If math.h has no long double support we can't rely 
   // on the math functions generating exponents outside 
   // the range of a double:
   //
   T min_val = (std::max)(
      tools::min_value(a), 
      static_cast<T>((std::numeric_limits<double>::min)()));
#else
   T min_val = tools::min_value(a);
#endif

   if((a != 0) && (b != 0))
   {
      // TODO: use isfinite:
      if(std::numeric_limits<T>::is_specialized && (b > (std::numeric_limits<T>::max)()))
      {
         if(a > (std::numeric_limits<T>::max)())
            return 0;
      }
      // if the result is denormalised, treat all denorms as equivalent:
      if(std::numeric_limits<T>::is_specialized)
      {
         if((a < min_val) && (a > 0))
            a = min_val;
         else if((a > -min_val) && (a < 0))
            a = -min_val;
         if((b < min_val) && (b > 0))
            b = min_val;
         else if((b > -min_val) && (b < 0))
            b = -min_val;
      }
      return (std::max)(fabs((a-b)/a), fabs((a-b)/b));
   }

   // handle special case where one or both are zero:
   if(!std::numeric_limits<T>::is_specialized)
      return fabs(a-b);
   if(fabs(a) < min_val)
      a = min_val;
   if(fabs(b) < min_val)
      b = min_val;
   return (std::max)(fabs((a-b)/a), fabs((a-b)/b));
}

template <class Seq>
void print_row(const Seq& row)
{
   for(unsigned i = 0; i < row.size(); ++i)
   {
      if(i)
         std::cout << ", ";
      std::cout << row[i];
   }
   std::cout << std::endl;
}

//
// Function test accepts an matrix of input values (probably a 2D boost::array)
// and calls two functors for each row in the array - one calculates a value
// to test, and one extracts the expected value from the array (or possibly
// calculates it at high precision).  The two functors are usually simple lambda
// expressions.
//
template <class A, class F1, class F2>
test_result<typename calculate_result_type<A>::value_type> test(const A& a, F1 test_func, F2 expect_func)
{
   typedef typename A::value_type         row_type;
   typedef typename row_type::value_type  value_type;

   test_result<value_type> result;
   result.worst_case = 0;

   for(unsigned i = 0; i < a.size(); ++i)
   {
      const row_type& row = a[i];
      value_type point = test_func(row);
      value_type expected = expect_func(row);
      value_type err = relative_error(point, expected);
#ifdef BOOST_INSTRUMENT
      if(err != 0) 
      {
         std::cout << row[0] << " " << err;  
         if(std::numeric_limits<value_type>::is_specialized)
         {
            std::cout << " (" << err / std::numeric_limits<value_type>::epsilon() << "eps)";
         }
         std::cout << std::endl;
      }
#endif
      if(!boost::math::isfinite(point) && boost::math::isfinite(expected))
      {
         std::cout << "CAUTION: Found non-finite result, when a finite value was expected at entry " << i << "\n";
         std::cout << "Found: " << point << " Expected " << expected << " Error: " << err << std::endl;
         print_row(row);
      }
      if(err > 0.5)
      {
         std::cout << "CAUTION: Gross error found at entry " << i << ".\n";
         std::cout << "Found: " << point << " Expected " << expected << " Error: " << err << std::endl;
         print_row(row);
      }
      result.stat.add(err);
      if((result.stat.max)() == err)
         result.worst_case = i;
   }
   return result;
}

} } } // namespaces

#endif

