//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SF_DIGAMMA_HPP
#define BOOST_MATH_SF_DIGAMMA_HPP

#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/evaluation_type.hpp>
#include <boost/math/tools/error_handling.hpp>
#include <boost/math/constants/constants.hpp>

namespace boost{
namespace math{
namespace detail{
//
// Begin by defining the smallest value for which it is safe to
// use the asymptotic expansion for digamma:
//
inline unsigned digamma_large_lim(const mpl::int_<0>*)
{  return 20;  }

inline unsigned digamma_large_lim(const void*)
{  return 10;  }
//
// Implementations of the asymptotic expansion come next,
// the coefficients of the series have been evaluated
// in advance at high precision, and the series truncated
// at the first term that's too small to effect the result.
// Note that the series becomes divergent after a while
// so truncation is very important.
//
// This first one gives 34-digit precision for x >= 20:
//
template <class T>
T digamma_imp_large(T x, const mpl::int_<0>*)
{
   using namespace std; // ADL of std functions.
   static const T P[] = {
      0.083333333333333333333333333333333333333333333333333L,
      -0.0083333333333333333333333333333333333333333333333333L,
      0.003968253968253968253968253968253968253968253968254L,
      -0.0041666666666666666666666666666666666666666666666667L,
      0.0075757575757575757575757575757575757575757575757576L,
      -0.021092796092796092796092796092796092796092796092796L,
      0.083333333333333333333333333333333333333333333333333L,
      -0.44325980392156862745098039215686274509803921568627L,
      3.0539543302701197438039543302701197438039543302701L,
      -26.456212121212121212121212121212121212121212121212L,
      281.4601449275362318840579710144927536231884057971L,
      -3607.510546398046398046398046398046398046398046398L,
      54827.583333333333333333333333333333333333333333333L,
      -974936.82385057471264367816091954022988505747126437L,
      20052695.796688078946143462272494530559046688078946L,
      -472384867.72162990196078431372549019607843137254902L,
      12635724795.916666666666666666666666666666666666667L
   };
   x -= 1;
   T result = log(x);
   result += 1 / (2 * x);
   T z = 1 / (x*x);
   result -= z * tools::evaluate_polynomial(P, z);
   return result;
}
//
// 19-digit precision for x >= 10:
//
template <class T>
T digamma_imp_large(T x, const mpl::int_<64>*)
{
   using namespace std; // ADL of std functions.
   static const T P[] = {
      0.083333333333333333333333333333333333333333333333333L,
      -0.0083333333333333333333333333333333333333333333333333L,
      0.003968253968253968253968253968253968253968253968254L,
      -0.0041666666666666666666666666666666666666666666666667L,
      0.0075757575757575757575757575757575757575757575757576L,
      -0.021092796092796092796092796092796092796092796092796L,
      0.083333333333333333333333333333333333333333333333333L,
      -0.44325980392156862745098039215686274509803921568627L,
      3.0539543302701197438039543302701197438039543302701L,
      -26.456212121212121212121212121212121212121212121212L,
      281.4601449275362318840579710144927536231884057971L,
   };
   x -= 1;
   T result = log(x);
   result += 1 / (2 * x);
   T z = 1 / (x*x);
   result -= z * tools::evaluate_polynomial(P, z);
   return result;
}
//
// 17-digit precision for x >= 10:
//
template <class T>
T digamma_imp_large(T x, const mpl::int_<53>*)
{
   using namespace std; // ADL of std functions.
   static const T P[] = {
      0.083333333333333333333333333333333333333333333333333L,
      -0.0083333333333333333333333333333333333333333333333333L,
      0.003968253968253968253968253968253968253968253968254L,
      -0.0041666666666666666666666666666666666666666666666667L,
      0.0075757575757575757575757575757575757575757575757576L,
      -0.021092796092796092796092796092796092796092796092796L,
      0.083333333333333333333333333333333333333333333333333L,
      -0.44325980392156862745098039215686274509803921568627L
   };
   x -= 1;
   T result = log(x);
   result += 1 / (2 * x);
   T z = 1 / (x*x);
   result -= z * tools::evaluate_polynomial(P, z);
   return result;
}
//
// 9-digit precision for x >= 10:
//
template <class T>
T digamma_imp_large(T x, const mpl::int_<24>*)
{
   using namespace std; // ADL of std functions.
   static const T P[] = {
      0.083333333333333333333333333333333333333333333333333L,
      -0.0083333333333333333333333333333333333333333333333333L,
      0.003968253968253968253968253968253968253968253968254L
   };
   x -= 1;
   T result = log(x);
   result += 1 / (2 * x);
   T z = 1 / (x*x);
   result -= z * tools::evaluate_polynomial(P, z);
   return result;
}
//
// Now follow rational approximations over the range [1,2].
//
// 35-digit precision:
//
template <class T>
T digamma_imp_1_2(T x, const mpl::int_<0>*)
{
   //
   // Now the approximation, we use the form:
   //
   // digamma(x) = (x - root) * (Y + R(x-1))
   //
   // Where root is the location of the positive root of digamma,
   // Y is a constant, and R is optimised for low absolute error
   // compared to Y.
   //
   // Max error found at 128-bit long double precision:  5.631e-35
   // Maximum Deviation Found (approximation error):     1.960e-35
   //
   static const float Y = 0.99558162689208984375F;

   static const T root1 = 1569415565.0 / 1073741824uL;
   static const T root2 = (381566830.0 / 1073741824uL) / 1073741824uL;
   static const T root3 = ((111616537.0 / 1073741824uL) / 1073741824uL) / 1073741824uL;
   static const T root4 = (((503992063.0 / 1073741824uL) / 1073741824uL) / 1073741824uL) / 1073741824uL;
   static const T root5 = 0.75231638452626400509999138382223723380394592598054e-36L;

   static const T P[] = {    
      0.25479851061131551526977464225335883773L,
      -0.18680627276785631913044260934720955868L,
      -0.80363716137044194938509411262041851893L,
      -0.67238768156656204735119811446899962752L,
      -0.26578204405928542571891729794823602027L,
      -0.0577913339744471190265540061417732562L,
      -0.0071502599302302344168886861743721325075L,
      -0.00048820246337241825839308387373932012069L,
      -0.16502375271306657440386132787211370224e-4L,
      -0.20455790127365005976211054106294513345e-6L,
      -0.92211014694806710961336881668647739535e-11L,
      0.48185872432754633773293738021874799566e-13L
   };
   static const T Q[] = {    
      1,
      2.6212362318426769040294336710114361821L,
      2.6854465068215025330325911346470503289L,
      1.4324621846378328535720514425511165882L,
      0.44127839367218238858671258622208841262L,
      0.081441923652895712841114777288711746235L,
      0.0089575839091639588314918592219010970046L,
      0.0005595920067347432979501367393705505528L,
      0.17654934314236241423226008246523175037e-4L,
      0.20719702289565546390545985045173039544e-6L
   };
   T g = x - root1;
   g -= root2;
   g -= root3;
   g -= root4;
   g -= root5;
   T r = tools::evaluate_polynomial(P, x-1) / tools::evaluate_polynomial(Q, x-1);
   T result = g * Y + g * r;

   return result;
}
//
// 19-digit precision:
//
template <class T>
T digamma_imp_1_2(T x, const mpl::int_<64>*)
{
   //
   // Now the approximation, we use the form:
   //
   // digamma(x) = (x - root) * (Y + R(x-1))
   //
   // Where root is the location of the positive root of digamma,
   // Y is a constant, and R is optimised for low absolute error
   // compared to Y.
   //
   // Max error found at 80-bit long double precision:   5.016e-20
   // Maximum Deviation Found (approximation error):     3.575e-20
   //
   static const float Y = 0.99558162689208984375F;

   static const T root1 = 1569415565.0 / 1073741824uL;
   static const T root2 = (381566830.0 / 1073741824uL) / 1073741824uL;
   static const T root3 = 0.9016312093258695918615325266959189453125e-19L;

   static const T P[] = {    
      0.254798510611315515235L,
      -0.314628554532916496608L,
      -0.665836341559876230295L,
      -0.314767657147375752913L,
      -0.0541156266153505273939L,
      -0.00289268368333918761452L
   };
   static const T Q[] = {    
      1,
      2.1195759927055347547L,
      1.54350554664961128724L,
      0.486986018231042975162L,
      0.0660481487173569812846L,
      0.00298999662592323990972L,
      -0.165079794012604905639e-5L,
      0.317940243105952177571e-7L
   };
   T g = x - root1;
   g -= root2;
   g -= root3;
   T r = tools::evaluate_polynomial(P, x-1) / tools::evaluate_polynomial(Q, x-1);
   T result = g * Y + g * r;

   return result;
}
//
// 18-digit precision:
//
template <class T>
T digamma_imp_1_2(T x, const mpl::int_<53>*)
{
   //
   // Now the approximation, we use the form:
   //
   // digamma(x) = (x - root) * (Y + R(x-1))
   //
   // Where root is the location of the positive root of digamma,
   // Y is a constant, and R is optimised for low absolute error
   // compared to Y.
   //
   // Maximum Deviation Found:               1.466e-18
   // At double precision, max error found:  2.452e-17
   //
   static const float Y = 0.99558162689208984F;

   static const T root1 = 1569415565.0 / 1073741824uL;
   static const T root2 = (381566830.0 / 1073741824uL) / 1073741824uL;
   static const T root3 = 0.9016312093258695918615325266959189453125e-19L;

   static const T P[] = {    
      0.25479851061131551L,
      -0.32555031186804491L,
      -0.65031853770896507L,
      -0.28919126444774784L,
      -0.045251321448739056L,
      -0.0020713321167745952L
   };
   static const T Q[] = {    
      1L,
      2.0767117023730469L,
      1.4606242909763515L,
      0.43593529692665969L,
      0.054151797245674225L,
      0.0021284987017821144L,
      -0.55789841321675513e-6L
   };
   T g = x - root1;
   g -= root2;
   g -= root3;
   T r = tools::evaluate_polynomial(P, x-1) / tools::evaluate_polynomial(Q, x-1);
   T result = g * Y + g * r;

   return result;
}
//
// 9-digit precision:
//
template <class T>
T digamma_imp_1_2(T x, const mpl::int_<24>*)
{
   //
   // Now the approximation, we use the form:
   //
   // digamma(x) = (x - root) * (Y + R(x-1))
   //
   // Where root is the location of the positive root of digamma,
   // Y is a constant, and R is optimised for low absolute error
   // compared to Y.
   //
   // Maximum Deviation Found:              3.388e-010
   // At float precision, max error found:  2.008725e-008
   //
   static const float Y = 0.99558162689208984f;
   static const T root = 1532632.0f / 1048576;
   static const T root_minor = static_cast<T>(0.3700660185912626595423257213284682051735604e-6L);
   static const T P[] = {    
      0.25479851023250261e0,
      -0.44981331915268368e0,
      -0.43916936919946835e0,
      -0.61041765350579073e-1
   };
   static const T Q[] = {    
      0.1e1,
      0.15890202430554952e1,
      0.65341249856146947e0,
      0.63851690523355715e-1
   };
   T g = x - root;
   g -= root_minor;
   T r = tools::evaluate_polynomial(P, x-1) / tools::evaluate_polynomial(Q, x-1);
   T result = g * Y + g * r;

   return result;
}

template <class T, class Tag>
T digamma_imp(T x, const Tag* t)
{
   //
   // This handles reflection of negative arguments, and all our
   // error handling, then forwards to the T-specific approximation.
   //
   using namespace std; // ADL of std functions.

   T result = 0;
   //
   // Check for negative arguments and use reflection:
   //
   if(x < 0)
   {
      // Reflect:
      x = 1 - x;
      // Argument reduction for tan:
      T remainder = x - floor(x);
      // Shift to negative if > 0.5:
      if(remainder > 0.5)
      {
         remainder -= 1;
      }
      //
      // check for evaluation at a negative pole:
      //
      if(remainder == 0)
      {
         return tools::pole_error<T>(BOOST_CURRENT_FUNCTION, 0, (1-x));
      }
      result = constants::pi<T>() / tan(constants::pi<T>() * remainder);
   }
   //
   // If we're above the lower-limit for the
   // asymptotic expansion then use it:
   //
   if(x >= digamma_large_lim(t))
   {
      result += digamma_imp_large(x, t);
   }
   else
   {
      //
      // If x > 2 reduce to the interval [1,2]:
      //
      while(x > 2)
      {
         x -= 1;
         result += 1/x;
      }
      //
      // If x < 1 use recurrance to shift to > 1:
      //
      if(x < 1)
      {
         result = -1/x;
         x += 1;
      }
      result += digamma_imp_1_2(x, t);
   }
   return result;
}

} // namespace detail

template <class T>
T digamma(T x)
{
   typedef typename tools::evaluation<typename remove_cv<T>::type>::type value_type;
   typedef typename mpl::if_c<
      (std::numeric_limits<T>::digits == 0)
      ||
      (std::numeric_limits<T>::digits > 64),
      mpl::int_<0>,
      typename mpl::if_c<
         (std::numeric_limits<T>::digits < 25),
         mpl::int_<24>,
         typename mpl::if_c<
            (std::numeric_limits<T>::digits < 54),
            mpl::int_<53>,
            mpl::int_<64>
         >::type
      >::type
   >::type tag_type;

   return tools::checked_narrowing_cast<typename remove_cv<T>::type>(detail::digamma_imp(
      static_cast<value_type>(x),
      static_cast<const tag_type*>(0)), BOOST_CURRENT_FUNCTION);
}

} // namespace math
} // namespace boost
#endif
