//  Copyright (c) 2007 John Maddock
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_COS_PI_HPP
#define BOOST_MATH_COS_PI_HPP

#include <cmath>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/real_cast.hpp>
#include <boost/math/constants/constants.hpp>

namespace boost{ namespace math{

template <class T>
T cos_pi(T x)
{
   using namespace std; // ADL of std names
   // cos of pi*x:
   bool invert = false;
   if(x < 0.5)
      return cos(constants::pi<T>() * x);
   if(x < 1)
   {
      x = -x;
   }

   T rem = floor(x);
   if(tools::real_cast<int>(rem) & 1)
      invert = !invert;
   rem = x - rem;
   if(rem > 0.5f)
   {
      rem = 1 - rem;
      invert = !invert;
   }
   if(rem == 0.5f)
      return 0;
   
   rem = cos(constants::pi<T>() * rem);
   return invert ? -rem : rem;
}

template <class T, class Policy>
inline T cos_pi(T x, const Policy&)
{
   return cos_pi(x);
}

} // namespace math
} // namespace boost
#endif
