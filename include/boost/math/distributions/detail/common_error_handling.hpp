// Copyright John Maddock 2006.

// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DISTRIBUTIONS_COMMON_ERROR_HANDLING_HPP
#define BOOST_MATH_DISTRIBUTIONS_COMMON_ERROR_HANDLING_HPP

#include <boost/math/tools/error_handling.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost{ namespace math{ namespace detail{

template <class RealType>
inline bool check_probability(const char* function, RealType const& prob, RealType* result)
{
   if((prob < 0) || (prob > 1) || !(boost::math::isfinite)(prob))
   {
      *result = tools::domain_error<RealType>(
         function, 
         "Probability argument is %1%, but must be >= 0 and <= 1 !", prob);
      return false;
   }
   return true;
}

template <class RealType>
inline bool check_df(const char* function, RealType const& df, RealType* result)
{
   if((df <= 0) || !(boost::math::isfinite)(df))
   {
      *result = tools::domain_error<RealType>(
         function, 
         "Degrees of freedom argument is %1%, but must be > 0 !", df);
      return false;
   }
   return true;
}


} // namespace detail
} // namespace math
} // namespace boost

#endif // BOOST_MATH_DISTRIBUTIONS_COMMON_ERROR_HANDLING_HPP
