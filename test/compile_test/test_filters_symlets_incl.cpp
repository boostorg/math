//  Copyright Nick Thompson 2026.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// Basic sanity check that header <boost/math/filters/symlets.hpp>
// #includes all the files that it needs to.
//
#include <boost/math/filters/symlets.hpp>

void compile_and_link_test()
{
   using namespace boost::math::filters;

   constexpr std::array<double, 4> da = symlet_scaling_filter<double, 2>();
   static_assert(da[0] == 0.48296291314453414337487159986444868381695241950420227520117153815521160699);
   constexpr std::array<float, 16> fa = symlet_scaling_filter<float, 8>();
   static_assert(fa.size() == 16);
}
