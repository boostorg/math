//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define L22
#include "../tools/ntl_rr_lanczos.hpp"
#include "../tools/ntl_rr_digamma.hpp"
#include <boost/math/bindings/rr.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>

#include <cmath>


boost::math::ntl::RR f(const boost::math::ntl::RR& x, int variant)
{
   static const boost::math::ntl::RR tiny = boost::math::tools::min_value<boost::math::ntl::RR>() * 64;
   switch(variant)
   {
   case 0:
      return boost::math::erfc(x) * x / exp(-x * x);
   case 1:
      return boost::math::erf(x);
   case 2:
      {
         boost::math::ntl::RR x_ = x == 0 ? 1e-80 : x;
      return boost::math::erf(x_) / x_;
      }
   case 3:
      {
         boost::math::ntl::RR y(x);
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
         return boost::lexical_cast<boost::math::ntl::RR>("0.42278433509846713939348790991759756895784066406008") / 3;
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
         boost::math::ntl::RR r1 = boost::lexical_cast<boost::math::ntl::RR>("0.57721566490153286060651209008240243104215933593992");
         boost::math::ntl::RR r2 = boost::lexical_cast<boost::math::ntl::RR>("0.42278433509846713939348790991759756895784066406008");
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
         boost::math::ntl::RR r1 = boost::lexical_cast<boost::math::ntl::RR>("0.57721566490153286060651209008240243104215933593992");
         boost::math::ntl::RR r2 = boost::lexical_cast<boost::math::ntl::RR>("0.42278433509846713939348790991759756895784066406008");
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
         boost::math::ntl::RR y = x;
         if(y == 0)
            y = boost::math::tools::epsilon<boost::math::ntl::RR>() / 64;
         return boost::math::erf_inv(y) / (y * (y+10));
      }
   case 8:
      {
         // 
         // erfc_inv in range [0.25, 0.5]
         // Use an y-offset of 0.25, and range [0, 0.25]
         // abs error, auto y-offset.
         //
         boost::math::ntl::RR y = x;
         if(y == 0)
            y = boost::lexical_cast<boost::math::ntl::RR>("1e-5000");
         return sqrt(-2 * log(y)) / boost::math::erfc_inv(y);
      }
   case 9:
      {
         boost::math::ntl::RR x2 = x;
         if(x2 == 0)
            x2 = boost::lexical_cast<boost::math::ntl::RR>("1e-5000");
         boost::math::ntl::RR y = exp(-x2*x2); // sqrt(-log(x2)) - 5;
         return boost::math::erfc_inv(y) / x2;
      }
   case 10:
      {
         //
         // Digamma over the interval [1,2], set x-offset to 1
         // and optimise for absolute error over [0,1].
         //
         int current_precision = boost::math::ntl::RR::precision();
         if(current_precision < 1000)
            boost::math::ntl::RR::SetPrecision(1000);
         //
         // This value for the root of digamma is calculated using our
         // differentiated lanczos approximation.  It agrees with Cody
         // to ~ 25 digits and to Morris to 35 digits.  See:
         // TOMS ALGORITHM 708 (Didonato and Morris).
         // and Math. Comp. 27, 123-127 (1973) by Cody, Strecok and Thacher.
         //
         //boost::math::ntl::RR root = boost::lexical_cast<boost::math::ntl::RR>("1.4616321449683623412626595423257213234331845807102825466429633351908372838889871");
         //
         // Actually better to calculate the root on the fly, it appears to be more
         // accurate: convergence is easier with the 1000-bit value, the approximation
         // produced agrees with functions.mathworld.com values to 35 digits even quite
         // near the root.
         //
         static boost::math::tools::eps_tolerance<boost::math::ntl::RR> tol(1000);
         static boost::uintmax_t max_iter = 1000;
         static const boost::math::ntl::RR root = boost::math::tools::bracket_and_solve_root(&boost::math::digamma, boost::math::ntl::RR(1.4), boost::math::ntl::RR(1.5), true, tol, max_iter).first;

         boost::math::ntl::RR x2 = x;
         double lim = 1e-65;
         if(fabs(x2 - root) < lim)
         {
            //
            // This is a problem area:
            // x2-root suffers cancellation error, so does digamma.
            // That gets compounded again when Remez calculates the error
            // function.  This cludge seems to stop the worst of the problems:
            //
            static const boost::math::ntl::RR a = boost::math::digamma(root - lim) / -lim;
            static const boost::math::ntl::RR b = boost::math::digamma(root + lim) / lim;
            boost::math::ntl::RR fract = (x2 - root + lim) / (2*lim);
            boost::math::ntl::RR r = (1-fract) * a + fract * b;
            std::cout << "In root area: " << r;
            return r;
         }
         boost::math::ntl::RR result =  boost::math::digamma(x2) / (x2 - root);
         if(current_precision < 1000)
            boost::math::ntl::RR::SetPrecision(current_precision);
         return result;
      }
   case 11:
      // expm1:
      if(x == 0)
      {
         static boost::math::ntl::RR lim = 1e-80;
         static boost::math::ntl::RR a = boost::math::expm1(-lim);
         static boost::math::ntl::RR b = boost::math::expm1(lim);
         static boost::math::ntl::RR l = (b-a) / (2 * lim);
         return l;
      }
      return boost::math::expm1(x) / x;
   case 12:
      // demo, and test case:
      return exp(x);
   case 13:
      // K(k):
      {
         // x = k^2
         boost::math::ntl::RR k2 = x;
         if(k2 > boost::math::ntl::RR(1) - 1e-40)
            k2 = boost::math::ntl::RR(1) - 1e-40;
         /*if(k < 1e-40)
            k = 1e-40;*/
         boost::math::ntl::RR p2 = boost::math::constants::pi<boost::math::ntl::RR>() / 2;
         return (boost::math::ellint_1(sqrt(k2))) / (p2 - boost::math::log1p(-k2));
      }
   case 14:
      // K(k)
      {
         // x = 1 - k^2
         boost::math::ntl::RR mp = x;
         if(mp < 1e-20)
            mp = 1e-20;
         boost::math::ntl::RR k = sqrt(1 - mp);
         static const boost::math::ntl::RR l4 = log(boost::math::ntl::RR(4));
         boost::math::ntl::RR p2 = boost::math::constants::pi<boost::math::ntl::RR>() / 2;
         return boost::math::ellint_1(k) / (l4 - log(mp));
      }
   case 15:
      // E(k)
      {
         // x = 1-k^2
         boost::math::ntl::RR z = 1 - x * log(x);
         return boost::math::ellint_2(sqrt(1-x)) / z;
      }
   }
   return 0;
}

void show_extra(
   const boost::math::tools::polynomial<boost::math::ntl::RR>& n, 
   const boost::math::tools::polynomial<boost::math::ntl::RR>& d, 
   const boost::math::ntl::RR& x_offset, 
   const boost::math::ntl::RR& y_offset, 
   int variant)
{
   switch(variant)
   {
   default:
      // do nothing here...
      ;
   }
}

