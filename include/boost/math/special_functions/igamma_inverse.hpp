//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_IGAMMA_INVERSE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_IGAMMA_INVERSE_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/tools/roots.hpp>

namespace boost{ namespace math{ namespace detail{

template <class T>
T estimate_inverse_s(T p, T q)
{
   //
   // Computation of the Incomplete Gamma Function Ratios and their Inverse
   // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
   // ACM Transactions on Mathematical Software, Vol. 12, No. 4, 
   // December 1986, Pages 377-393.
   //
   // See equation 32.
   //
   using namespace std;
   T t;
   if(p < 0.5)
   {
      t = sqrt(-2 * log(p));
   }
   else
   {
      t = sqrt(-2 * log(q));
   }
   static const double a[4] = { 3.31125922108741, 11.6616720288968, 4.28342155967104, 0.213623493715853 };
   static const double b[5] = { 1, 6.61053765625462, 6.40691597760039, 1.27364489782223, 0.3611708101884203e-1 };
   T s = t - tools::evaluate_polynomial(a, t, 4) / tools::evaluate_polynomial(b, t, 5);
   if(p < 0.5)
      s = -s;
   return s;
}

template <class T>
T didonato_SN(T a, T x, unsigned N, T tolerance = 0)
{
   //
   // Computation of the Incomplete Gamma Function Ratios and their Inverse
   // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
   // ACM Transactions on Mathematical Software, Vol. 12, No. 4, 
   // December 1986, Pages 377-393.
   //
   // See equation 34.
   //
   T sum = 1;
   if(N >= 1)
   {
      T partial = x / (a + 1);
      sum += partial;
      for(unsigned i = 2; i <= N; ++i)
      {
         partial *= x / (a + i);
         sum += partial;
         if(partial < tolerance)
            break;
      }
   }
   return sum;
}

template <class T>
T didonato_FN(T p, T a, T x, unsigned N, T tolerance)
{
   //
   // Computation of the Incomplete Gamma Function Ratios and their Inverse
   // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
   // ACM Transactions on Mathematical Software, Vol. 12, No. 4, 
   // December 1986, Pages 377-393.
   //
   // See equation 34.
   //
   using namespace std;
   T u = log(p) + boost::math::lgamma(a + 1);
   return exp((u + x - log(didonato_SN(a, x, N, tolerance))) / a);
}

template <class T>
T estimate_inverse_gamma(T a, T p, T q)
{
   //
   // In order to understand what's going on here, you will
   // need to refer to:
   //
   // Computation of the Incomplete Gamma Function Ratios and their Inverse
   // ARMIDO R. DIDONATO and ALFRED H. MORRIS, JR.
   // ACM Transactions on Mathematical Software, Vol. 12, No. 4, 
   // December 1986, Pages 377-393.
   //
   using namespace std;

   T result;

   if(a == 1)
      result = -log(q);
   else if(a < 1)
   {
      T g = boost::math::tgamma(a);
      T b = q * g;
      if((b > 0.6) || ((b >= 0.45) && (a >= 0.3)))
      {
         // DiDonato & Morris Eq 21:
         //
         // There is a slight variation from DiDonato and Morris here:
         // the first form given here is unstable when p is close to 1,
         // making it impossible to compute the inverse of Q(a,x) for small
         // q.  Fortunately the second form works perfectly well in this case.
         //
         T u;
         if((b * q > 1e-8) && (q > 1e-5))
         {
            u = pow(p * g * a, 1 / a);
         }
         else
         {
            u = exp((-q / a) - constants::euler<T>());
         }
         result = u / (1 - (u / (a + 1)));
      }
      else if((a < 0.3) && (b >= 0.35))
      {
         // DiDonato & Morris Eq 22:
         T t = exp(-constants::euler<T>() - b);
         T u = t * exp(t);
         result = t * exp(u);
      }
      else if((b > 0.15) || (a >= 0.3))
      {
         // DiDonato & Morris Eq 23:
         T y = -log(b);
         T u = y - (1 - a) * log(y);
         result = y - (1 - a) * log(u) - log(1 + (1 - a) / (1 + u));
      }
      else if (b > 0.1)
      {
         // DiDonato & Morris Eq 24:
         T y = -log(b);
         T u = y - (1 - a) * log(y);
         result = y - (1 - a) * log(u) - log((u * u + 2 * (3 - a) * u + (2 - a) * (3 - a)) / (u * u + (5 - a) * u + 2));
      }
      else
      {
         // DiDonato & Morris Eq 25:
         T y = -log(b);
         T c1 = (a - 1) * log(y);
         T c1_2 = c1 * c1;
         T c1_3 = c1_2 * c1;
         T c1_4 = c1_2 * c1_2;
         T a_2 = a * a;
         T a_3 = a_2 * a;

         T c2 = (a - 1) * (1 + c1);
         T c3 = (a - 1) * (-(c1_2 / 2) + (a - 2) * c1 + (3 * a - 5) / 2);
         T c4 = (a - 1) * ((c1_3 / 3) - (3 * a - 5) * c1_2 / 2 + (a_2 - 6 * a + 7) * c1 + (11 * a_2 - 46 * a + 47) / 6);
         T c5 = (a - 1) * (-(c1_4 / 4) 
                           + (11 * a - 17) * c1_3 / 6 
                           + (-3 * a_2 + 13 * a -13) * c1_2 
                           + (2 * a_3 - 25 * a_2 + 72 * a - 61) * c1 / 2
                           + (25 * a_3 - 195 * a_2 + 477 * a - 379) / 12);

         T y_2 = y * y;
         T y_3 = y_2 * y;
         T y_4 = y_2 * y_2;
         result = y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4);
      }
   }
   else
   {
      // DiDonato and Morris Eq 31:
      T s = estimate_inverse_s(p, q);

      T s_2 = s * s;
      T s_3 = s_2 * s;
      T s_4 = s_2 * s_2;
      T s_5 = s_4 * s;
      T ra = sqrt(a);
      
      T w = a + s * ra + (s * s -1) / 3;
      w += (s_3 - 7 * s) / (36 * ra);
      w -= (3 * s_4 + 7 * s_2 - 16) / (810 * a);
      w += (9 * s_5 + 256 * s_3 - 433 * s) / (38880 * a * ra);

      if((a >= 500) && (fabs(1 - w / a) < 1e-6))
      {
         result = w;
      }
      else if (p > 0.5)
      {
         if(w < 3 * a)
         {
            result = w;
         }
         else
         {
            T D = (std::max)(T(2), a * (a - 1));
            T lg = boost::math::lgamma(a);
            T lb = log(q) + lg;
            if(lb < -D * 2.3)
            {
               // DiDonato and Morris Eq 25:
               T y = -lb;
               T c1 = (a - 1) * log(y);
               T c1_2 = c1 * c1;
               T c1_3 = c1_2 * c1;
               T c1_4 = c1_2 * c1_2;
               T a_2 = a * a;
               T a_3 = a_2 * a;

               T c2 = (a - 1) * (1 + c1);
               T c3 = (a - 1) * (-(c1_2 / 2) + (a - 2) * c1 + (3 * a - 5) / 2);
               T c4 = (a - 1) * ((c1_3 / 3) - (3 * a - 5) * c1_2 / 2 + (a_2 - 6 * a + 7) * c1 + (11 * a_2 - 46 * a + 47) / 6);
               T c5 = (a - 1) * (-(c1_4 / 4) 
                                 + (11 * a - 17) * c1_3 / 6 
                                 + (-3 * a_2 + 13 * a -13) * c1_2 
                                 + (2 * a_3 - 25 * a_2 + 72 * a - 61) * c1 / 2
                                 + (25 * a_3 - 195 * a_2 + 477 * a - 379) / 12);

               T y_2 = y * y;
               T y_3 = y_2 * y;
               T y_4 = y_2 * y_2;
               result = y + c1 + (c2 / y) + (c3 / y_2) + (c4 / y_3) + (c5 / y_4);
            }
            else
            {
               // DiDonato and Morris Eq 33:
               T u = -lb + (a - 1) * log(w) - log(1 + (1 - a) / (1 + w));
               result = -lb + (a - 1) * log(u) - log(1 + (1 - a) / (1 + u));
            }
         }
      }
      else
      {
         // DiDonato and Morris Eq 35:
         T z = didonato_FN(p, a, w, 0, T(0));
         z = didonato_FN(p, a, z, 2, T(0));
         z = didonato_FN(p, a, z, 2, T(0));
         z = didonato_FN(p, a, z, 3, T(0));

         if((z <= 0.01 * (a + 1)) || (z > 0.7 * (a + 1)))
         {
            result = z;
         }
         else
         {
            // DiDonato and Morris Eq 36:
            T zb = didonato_FN(p, a, z, 100, T(1e-4));
            T u = log(p) + boost::math::lgamma(a + 1);
            result = zb * (1 - (a * log(zb) - zb - u + log(didonato_SN(a, z, 100, T(1e-4)))) / (a - zb));
         }
      }
   }
   return result;
}

template <class T>
struct gamma_P_inverse_func
{
   gamma_P_inverse_func(T a_, T p_, bool inv) : a(a_), p(p_), invert(inv) 
   {
      //
      // If p is too near 1 then P(x) - p suffers from cancellation
      // errors causing our root-finding algorithms to "thrash", better
      // to invert in this case and calculate Q(x) - (1-p) instead.
      //
      // Of course if p is *very* close to 1, then the answer we get will
      // be inaccurate anyway (because there's not enough information in p)
      // but at least we will converge on the (inaccurate) answer quickly.
      //
      if(p > 0.9)
      {
         p = 1 - p;
         invert = !invert;
      }
   }

   std::tr1::tuple<T, T, T> operator()(const T& x)const
   {
      //
      // Calculate P(x) - p and the first two derivates, or if the invert
      // flag is set, then Q(x) - q and it's derivatives.
      //
      typedef typename lanczos::lanczos_traits<T>::value_type value_type;
      typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;

      T f = !invert ? boost::math::gamma_P(a, x) : boost::math::gamma_Q(a, x);
      T f1 = static_cast<T>(boost::math::detail::regularised_gamma_prefix(value_type(a), value_type(x), evaluation_type()));
      T f2;
      if((x < 1) && (tools::max_value(f1) * x < fabs(f1)))
      {
         // overflow:
         f1 = tools::max_value(f1) / 2;
         f2 = -f1;
      }
      else
      {
         f1 /= x;
         T div = (a - x - 1) / x;
         f2 = f1;
         if((fabs(div) > 1) && (tools::max_value(f2) / fabs(div) < f2))
         {
            // overflow:
            f2 = tools::max_value(f1) / 2;
         }
         else
         {
            f2 *= div;
         }
      }
      if(invert)
      {
         f1 = -f1;
         f2 = -f2;
      }

      return std::tr1::make_tuple(f - p, f1, f2);
   }
private:
   T a, p;
   bool invert;
};

} // namespace detail

template <class T>
T gamma_P_inv(T a, T p)
{
   if(p == 1)
      return tools::max_value(p);
   if(p == 0)
      return 0;
   T guess = detail::estimate_inverse_gamma(a, p, 1 - p);
   if((guess == 0) && (p != 0))
      guess = tools::min_value(guess);
   return tools::halley_iterate(
      detail::gamma_P_inverse_func<T>(a, p, false), 
      guess,
      T(0),
      tools::max_value(a),
      (tools::digits(a) * 2) / 3);
}

template <class T>
T gamma_Q_inv(T a, T q)
{
   if(q == 0)
      return tools::max_value(q);
   if(q == 1)
      return 0;
   T guess = detail::estimate_inverse_gamma(a, 1 - q, q);
   if((guess == 0) && (q != 1))
      guess = tools::min_value(guess);
   return tools::halley_iterate(
      detail::gamma_P_inverse_func<T>(a, q, true), 
      guess,
      T(0),
      tools::max_value(a),
      (tools::digits(a) * 2) / 3);
}

}} // namespaces

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_IGAMMA_INVERSE_HPP

