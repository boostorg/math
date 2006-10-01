//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define L22
#include "../tools/ntl_rr_lanczos.hpp"
#include "../tools/ntl_rr_digamma.hpp"
#include <boost/math/tools/ntl.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>

#include <cmath>


NTL::RR f(const NTL::RR& x, int variant)
{
   static const NTL::RR tiny = boost::math::tools::min_value<NTL::RR>() * 64;
   switch(variant)
   {
   case 0:
      return boost::math::expm1(x);
   case 1:
      return boost::math::log1p(x) - x;
   case 2:
      return boost::math::erf(x) / x - 1.125;
   case 3:
      {
         NTL::RR y(x);
         if(y == 0) 
            y += tiny;
         return boost::math::lgamma(y+2) / y - 0.5;
      }
   case 4:
      //
      // lgamma in the range [2,3], use:
      //
      // lgamma(x) = (x-2) * (x + 1) * (c + R(x - 2))
      //
      // Works well at 80-bit long double precision, but doesn't
      // stretch to 128-bit precision.
      //
      if(x == 0)
      {
         return boost::lexical_cast<NTL::RR>("0.42278433509846713939348790991759756895784066406008") / 3;
      }
      return boost::math::lgamma(x+2) / (x * (x+3));
   case 5:
      {
         //
         // lgamma in the range [1,2], use:
         //
         // lgamma(x) = (x - 1) * (x - 2) * (c + R(x - 1))
         //
         // works well over [1, 1.5] but not near 2 :-(
         //
         NTL::RR r1 = boost::lexical_cast<NTL::RR>("0.57721566490153286060651209008240243104215933593992");
         NTL::RR r2 = boost::lexical_cast<NTL::RR>("0.42278433509846713939348790991759756895784066406008");
         if(x == 0)
         {
            return r1;
         }
         if(x == 1)
         {
            return r2;
         }
         return boost::math::lgamma(x+1) / (x * (x - 1));
      }
   case 6:
      {
         //
         // lgamma in the range [1.5,2], use:
         //
         // lgamma(x) = (2 - x) * (1 - x) * (c + R(2 - x))
         //
         // works well over [1.5, 2] but not near 1 :-(
         //
         NTL::RR r1 = boost::lexical_cast<NTL::RR>("0.57721566490153286060651209008240243104215933593992");
         NTL::RR r2 = boost::lexical_cast<NTL::RR>("0.42278433509846713939348790991759756895784066406008");
         if(x == 0)
         {
            return r2;
         }
         if(x == 1)
         {
            return r1;
         }
         return boost::math::lgamma(2-x) / (x * (x - 1));
      }
   case 7:
      {
         //
         // erf_inv in range [0, 0.5]
         //
         NTL::RR y = x;
         if(y == 0)
            y = boost::math::tools::epsilon<NTL::RR>() / 64;
         return boost::math::erf_inv(y) / (y * (y+10));
      }
   case 8:
      {
         // 
         // erfc_inv in range [0.25, 0.5]
         // Use an y-offset of 0.25, and range [0, 0.25]
         // abs error, auto y-offset.
         //
         NTL::RR y = x;
         if(y == 0)
            y = boost::lexical_cast<NTL::RR>("1e-5000");
         return sqrt(-2 * log(y)) / boost::math::erfc_inv(y);
      }
   case 9:
      {
         NTL::RR x2 = x;
         if(x2 == 0)
            x2 = boost::lexical_cast<NTL::RR>("1e-5000");
         NTL::RR y = exp(-x2*x2); // sqrt(-log(x2)) - 5;
         return boost::math::erfc_inv(y) / x2;
      }
   case 10:
      {
         int current_precision = NTL::RR::precision();
         if(current_precision < 1000)
            NTL::RR::SetPrecision(1000);
         //
         // This value for the root of digamma is calculated using our
         // differentiated lanczos approximation.  It agrees with Cody
         // to ~ 25 digits and to Morris to 35 digits.  See:
         // TOLMS ALGORITHM 708 (Didonato and Morris).
         // and Math. Comp. 27, 123-127 (1973) by Cody, Strecok and Thacher.
         //
         NTL::RR root = boost::lexical_cast<NTL::RR>("1.4616321449683623412626595423257213234331845807102825466429633351908372838889871");

         NTL::RR x2 = x;
         if(x2 == root)
            x2 += 1e-40;
         NTL::RR result =  boost::math::digamma(x2) / (x2 - root);
         if(current_precision < 1000)
            NTL::RR::SetPrecision(current_precision);
         return result;
      }
   }
   return 0;
}

void show_extra(
   const boost::math::tools::polynomial<NTL::RR>& n, 
   const boost::math::tools::polynomial<NTL::RR>& d, 
   const NTL::RR& x_offset, 
   const NTL::RR& y_offset, 
   int variant)
{
   switch(variant)
   {
   default:
      // do nothing here...
      ;
   }
}

