

#include "../tools/ntl_rr_lanczos.hpp"
#include <boost/math/tools/ntl.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/erf.hpp>

#include <cmath>

template <class T>
T f0(T x, int variant)
{
   static const T tiny = boost::math::tools::min_value<T>() * 64;
   switch(variant)
   {
   case 0:
      return boost::math::expm1(x);
   case 1:
      return boost::math::log1p(x) - x;
   case 2:
      return boost::math::erf(x) / x - 1.125;
   case 3:
      if(x == 0) 
         x += tiny;
      return boost::math::lgamma(x+2) / x - 0.5;
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
         return boost::lexical_cast<T>("0.42278433509846713939348790991759756895784066406008") / 3;
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
         NTL::RR r1 = boost::lexical_cast<T>("0.57721566490153286060651209008240243104215933593992");
         NTL::RR r2 = boost::lexical_cast<T>("0.42278433509846713939348790991759756895784066406008");
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
         NTL::RR r1 = boost::lexical_cast<T>("0.57721566490153286060651209008240243104215933593992");
         NTL::RR r2 = boost::lexical_cast<T>("0.42278433509846713939348790991759756895784066406008");
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
      //
      // erf_inv in range [0, 0.5]
      //
      if(x == 0)
         x = boost::math::tools::epsilon<T>() / 64;
      return boost::math::erf_inv(x) / x;
   case 8:
      // 
      //erfc_inv in range [0.5, 1]
      //
      if(x == 0)
         return x = boost::lexical_cast<T>("1e-5000");
      return sqrt(-2 * log(x)) / boost::math::erfc_inv(x);
   case 9:
      if(x == 0)
         return 0;
      return 1 / boost::math::erfc_inv(x);
   }
   return 0;
}


NTL::RR f(const NTL::RR& x, int variant)
{
   return f0<NTL::RR>(x, variant);
}

void show_extra(
   const boost::math::tools::polynomial<NTL::RR>& n, 
   const boost::math::tools::polynomial<NTL::RR>& d, 
   const NTL::RR& x_offset, 
   const NTL::RR& y_offset, 
   int variant)
{
}

