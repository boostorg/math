//  Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_BINOMIAL_HPP
#define BOOST_MATH_SF_BINOMIAL_HPP

#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/tools/error_handling.hpp>

namespace boost{ namespace math{

template <class T>
T binomial_coefficient(unsigned n, unsigned k)
{
   using namespace std;
   if(k > n)
      return tools::domain_error<T>(
         BOOST_CURRENT_FUNCTION, 
         "The binomial coefficient is undefined for k > n, but got k = %1%.",
         k);
   T result;
   if((k == 0) || (k == n))
      return 1;
   if((k == 1) || (k == n-1))
      return n;

   if(n <= max_factorial<T>::value)
   {
      // Use fast table lookup:
      result = unchecked_factorial<T>(n);
      result /= unchecked_factorial<T>(n-k);
      result /= unchecked_factorial<T>(k);
   }
   else
   {
      // Use the beta function:
      if(k < n - k)
         result = k * beta(static_cast<T>(k), static_cast<T>(n-k+1));
      else
         result = (n - k) * beta(static_cast<T>(k+1), static_cast<T>(n-k));
      if(result == 0)
         return tools::overflow_error<T>(BOOST_CURRENT_FUNCTION, 0);
      result = 1 / result;
   }
   // convert to nearest integer:
   return ceil(result - 0.5f);
}
//
// Type float can only store the first 35 factorials, in order to
// increase the chance that we can use a table driven implementation
// we'll promote to double:
//
template <>
inline float binomial_coefficient<float>(unsigned n, unsigned k)
{
   return tools::checked_narrowing_cast<float>(binomial_coefficient<double>(n, k), BOOST_CURRENT_FUNCTION);
}

} // namespace math
} // namespace boost


#endif // BOOST_MATH_SF_BINOMIAL_HPP

