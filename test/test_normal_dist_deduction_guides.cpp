//  (C) Copyright John Maddock 2022.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/math/distributions/normal.hpp>

int main()
{
   using namespace boost::math;

   normal_distribution n;
   static_assert(std::is_same<decltype(n)::value_type, double>::value);
   normal_distribution n2(2);
   static_assert(std::is_same<decltype(n2)::value_type, double>::value);
   normal_distribution n3(2, 3);
   static_assert(std::is_same<decltype(n3)::value_type, double>::value);
}
