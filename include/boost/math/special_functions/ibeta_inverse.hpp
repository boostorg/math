//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_IBETA_INVERSE_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_IBETA_INVERSE_HPP

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/roots.hpp>

namespace boost{ namespace math{ namespace detail{

template <class T>
struct temme_root_finder
{
   temme_root_finder(const T t_, const T a_) : t(t_), a(a_) {}

   std::tr1::tuple<T, T> operator()(const T& x)
   {
      T p = pow(1-x, a);
      T f = x * p - t;
      T f1 = p - x * a * a * p / (1 - x);
      return std::tr1::make_tuple(f, f1);
   }
private:
   T t, a;
};

template <class T>
struct temme_root_finder2
{
   temme_root_finder2(const T t_, const T a_) : t(t_), a(a_) {}

   std::tr1::tuple<T, T> operator()(const T& x)
   {
      T f = log(x) + a * log(1 - x) + t;
      T f1 = (1 / x) - (a / (1 - x));
      return std::tr1::make_tuple(f, f1);
   }
private:
   T t, a;
};

template <class T>
T estimate_inverse_beta(T a, T b, T z)
{
   using namespace std;
   if(z == 1)
      return 1;
   if(z == 0)
      return 0;

   if(a + b > 5)
   {
      if((fabs(a-b) < 5) && ((std::min)(a, b) > 5))
      {
         const T r2 = sqrt(T(2));
         bool invert = false;
         if(a < b)
         {
            std::swap(a, b);
            invert = true;
         }
         //
         // get the first approximation for nu from the inverse
         // error function:
         //
         T nu0 = boost::math::erfc_inv(2 * z);
         nu0 /= -sqrt(a / 2);
         //
         // calculate powers:
         //
         T B = a - b;
         T nu0_2 = nu0 * nu0;
         T nu0_3 = nu0_2 * nu0;
         T nu0_4 = nu0_2 * nu0_2;
         T nu0_5 = nu0_3 * nu0_3;
         T nu0_6 = nu0_3 * nu0_3;
         T B_2 = B * B;
         T B_3 = B_2 * B;
         //
         // Calculate correction terms:
         //
         T e1 = -B * r2 / 2 + (1 - 2 * B) * nu0 / 8;
         e1 -= B * r2 * nu0_2 / 48;
         e1 -= nu0_3 / 192;
         e1 -= B * r2 * nu0_4 / 3840;

         T e2 = B * r2 * (3 * B - 2) / 12;
         e2 += (20 * B_2 - 12 * B + 1) * nu0 / 128;
         e2 += B * r2 * (20 * B - 1) * nu0_2 / 960;
         e2 += (16 * B_2 + 30 * B - 15) * nu0_3 / 4608;
         e2 += B * r2 * (21 * B + 32) * nu0_4 / 53760;
         e2 += (-32 * B_2 + 63) * nu0_5 / 368640;
         e2 -= B * r2 * (120 * B + 17) * nu0_6 / 25804480;

         T e3 = B * r2 * (-75 * B_2 + 80 * B - 16) / 480;
         e3 += (-1080 * B_3 + 868 * B_2 - 90 * B - 45) * nu0 / 9216;
         e3 += B * r2 * (-1190 * B_2 + 84 * B + 373) * nu0_2 / 53760;
         e3 += (-2240 * B_3 - 2508 * B_2 + 2100 * B - 165) * nu0_3 / 368640;
         //
         // Bring them together to get a final estimate for nu:
         //
         T nu = nu0 + e1 / a + e2 / (a * a) + e3 / (a * a * a);
         //
         // now we need to convert nu to x, by solving the appropriate
         // quadratic equation:
         //
         T nu_2 = nu * nu;
         T c = -exp(-nu_2 / 2);
         T x = (1 + nu * sqrt((1 + c) / nu_2)) / 2;
         BOOST_ASSERT(x >= 0);
         BOOST_ASSERT(x <= 1);
#ifdef BOOST_INSTRUMENT
         std::cout << "Estimating x with Temme method 1: " << x << std::endl;
#endif
         return x;
      }
      T r = a + b;
      T theta = asin(sqrt(a / r));
      if((theta > 0.7) && (theta < 0.8))
      {
         T nu0 = boost::math::erfc_inv(2 * z);
         nu0 /= -sqrt(r / 2);

         T s = sin(theta);
         T c = cos(theta);
         T s_2 = s * s;
         T s_3 = s_2 * s;
         T s_4 = s_2 * s_2;
         T s_5 = s_3 * s_2;
         T s_6 = s_3 * s_3;
         T s_7 = s_4 * s_3;
         T s_8 = s_4 * s_4;
         T s_9 = s_4 * s_5;
         T s_10 = s_5 * s_5;
         T s_12 = s_6 * s_6;
         T s_14 = s_7 * s_7;
         T c_2 = c * c;
         T c_3 = s_2 * s;
         T c_4 = c_2 * c_2;
         T c_5 = c_2 * c_3;
         T c_6 = c_3 * c_3;
         T c_7 = c_4 * c_3;
         T nu_2 = nu0 * nu0;
         T nu_3 = nu_2 * nu0;
         T nu_4 = nu_2 * nu_2;

         T e1 = (2 * s_2 - 1) / (3 * s * c);
         e1 -= (5 * s_4 - 5 * s_2 - 1) * nu0 / (36 * s_2 * c_2);
         e1 += (46 * s_6 - 69 * s_4 + 21 * s_2 + 1) * nu_2 / (1620 * s_3 * c_3);
         e1 -= (-2 * s_2 - 62 * s_6 + 31 * s_8 + 33 * s_4 + 7) * nu_3 / (6480 * s_4 * c_4);
         e1 += (88 * s_6 - 52 * s_2 - 115 * s_8 + 46 * s_10 - 17 * s_4 + 25) * nu_4 / (90720 * s_5 * c_5);

         T e2 = -(52 * s_6 - 78 * s_4 + 12 * s_2 + 7) / (405 * s_3 * c_3);
         e2 += (2 * s_2 - 370 * s_6 + 185 * s_8 + 183 * s_4 - 7) * nu0 / (2592 * s_4 * c_4);
         e2 -= (776 * s_2 + 10240 * s_6 - 13525 * s_8 - 533 + 5410 * s_10 - 1835 * s_4) * nu_2 / (204120 * s_5 * c_5);
         e2 += (3747 * s_2 + 15071 * s_12 - 15821 * s_6 + 45588 * s_8 - 45213 * s_10 - 3372 * s_4 - 1579) * nu_3 / (2099520 * s_6 * c_6);

         T e3 = (3704 * s_10 - 9260 * s_8 + 6686 * s_6 - 769 * s_4 - 1259 * s_2 + 449) / (102060 * s_5 * c_5);
         e3 -= (750479 * s_12 - 151557 * s_2 - 727469 * s_6 + 2239932 * s_8 - 2251437 * s_10 + 140052 * s_4 + 63149) * nu0 / (20995200 * s_6 * c_6);
         e3 += (729754 * s_14 - 78755 * s_2 - 2554139 * s_12 + 146879 * s_6 - 1602610 * s_8 + 3195183 * s_10 + 105222 * s_4 + 29233) * nu_2 / (36741600 * s_7 * c_7);

         T nu = nu0 + e1 / r + e2 / (r*r) + e3 / (r*r*r);
         //
         // Now back solve for x:
         //
         T x;
         if(fabs(nu) < 0.7)
         {
            nu_2 = nu * nu;
            nu_3 = nu_2 * nu;
            nu_4 = nu_2 * nu_2;
            x = s_2 + s * c * nu + (1 - 2 * s_2) * nu_2 / 3 + (13 * s_4 - 13 * s_2 + 1) * nu_3 / (36 * s * c) + (46 * s_6 - 69 * s_4 + 21 * s_2 + 1) * nu_4 / (270 * s_2 * c_2);
#ifdef BOOST_INSTRUMENT
            std::cout << "Estimating x with Temme method 2 (small nu): " << x << std::endl;
#endif
         }
         else
         {
            T alpha = c_2 / s_2;
            T u = exp((-(nu * nu) / (2 * s_2) + log(s_2) + c_2 * log(c_2) / s_2));
            T u_2 = u * u;
            T u_3 = u_2 * u;
            T u_4 = u_2 * u_2;
            T u_5 = u_3 * u_2;

            x = u + alpha * u + 3 * alpha * (3 * alpha + 1) * u_3 / 6;
            x += 4 * alpha * (4 * alpha + 1) * (4 * alpha + 2) * u_4 / 24;
            x += 5 * alpha * (5 * alpha + 1) * (5 * alpha + 2) * (5 * alpha + 3) * u_5 / 120;
            //
            // This approximation to x so far is straight out of Temme's paper
            // but isn't all that accurate, need to clean it up with some
            // Newton-Raphson iterations to stand any chance of being close
            // to the incomplete Beta inverse.  This is still *much* cheaper
            // than iterating on the full incomplete beta:
            //
            x = tools::newton_raphson_iterate(
               temme_root_finder<T>(u, alpha), x, T(0), T(1), (2 * tools::digits(x)) / 3);
            if((x - s_2) * nu < 0)
               x = 1 - x;
#ifdef BOOST_INSTRUMENT
            std::cout << "Estimating x with Temme method 2 (large nu): " << x << std::endl;
#endif
         }
         return x;
      }
      /*
      if(b > a)
         return 1 - estimate_inverse_beta(b, a, 1 - z); */
      if(1 /*a > b*/)
      {
         T nu0;
         if(z > 0.5)
            nu0 = boost::math::gamma_P_inv(b, 1 - z);
         else
            nu0 = boost::math::gamma_Q_inv(b, z);
         nu0 /= a;
         T mu = b / a;
         T w = sqrt(1 + mu);
         T w_2 = w * w;
         T w_3 = w_2 * w;
         T w_4 = w_2 * w_2;
         T w_5 = w_3 * w_2;
         T w_6 = w_3 * w_3;
         T w_7 = w_4 * w_3;
         T w_8 = w_4 * w_4;
         T w_9 = w_5 * w_4;
         T w_10 = w_5 * w_5;
         T d = nu0 - mu;
         T d_2 = d * d;
         T d_3 = d_2 * d;
         T d_4 = d_2 * d_2;
         T w1 = w + 1;
         T w1_2 = w1 * w1;
         T w1_3 = w1 * w1_2;
         T w1_4 = w1_2 * w1_2;

         T e1 = (w + 2) * (w - 1) / (3 * w);
         e1 += (w_3 + 9 * w_2 + 21 * w + 5) * d / (36 * w_2 * w1);
         e1 -= (w_4 - 13 * w_3 + 69 * w_2 + 167 * w + 46) * d_2 / (1620 * w1_2 * w_3);
         e1 -= (7 * w_5 + 21 * w_4 + 70 * w_3 + 26 * w_2 - 93 * w - 31) * d_3 / (6480 * w1_3 * w_4);
         e1 -= (75 * w_6 + 202 * w_5 + 188 * w_4 - 888 * w_3 - 1345 * w_2 + 118 * w + 138) * d_4 / (272160 * w1_4 * w_5);

         T e2 = (28 * w_4 + 131 * w_3 + 402 * w_2 + 581 * w + 208) * (w - 1) / (1620 * w1 * w_3);
         e2 -= (35 * w_6 - 154 * w_5 - 623 * w_4 - 1636 * w_3 - 3983 * w_2 - 3514 * w - 925) * d / (12960 * w1_2 * w_4);
         e2 -= (2132 * w_7 + 7915 * w_6 + 16821 * w_5 + 35066 * w_4 + 87490 * w_3 + 141183 * w_2 + 95993 * w + 21640) * d_2  / (816480 * w_5 * w1_3);
         e2 -= (11053 * w_8 + 53308 * w_7 + 117010 * w_6 + 163924 * w_5 + 116188 * w_4 - 258428 * w_3 - 677042 * w_2 - 481940 * w - 105497) * d_3 / (14696640 * w1_4 * w_6);

         T e3 = -((3592 * w_7 + 8375 * w_6 - 1323 * w_5 - 29198 * w_4 - 89578 * w_3 - 154413 * w_2 - 116063 * w - 29632) * (w - 1)) / (816480 * w_5 * w1_2);
         e3 -= (442043 * w_9 + 2054169 * w_8 + 3803094 * w_7 + 3470754 * w_6 + 2141568 * w_5 - 2393568 * w_4 - 19904934 * w_3 - 34714674 * w_2 - 23128299 * w - 5253353) * d / (146966400 * w_6 * w1_3);
         e3 -= (116932 * w_10 + 819281 * w_9 + 2378172 * w_8 + 4341330 * w_7 + 6806004 * w_6 + 10622748 * w_5 + 18739500 * w_4 + 30651894 * w_3 + 30869976 * w_2 + 15431867 * w + 2919016) * d_2 / (146966400 * w1_4 * w_7);

         T nu = nu0 + e1 / a + e2 / (a * a) + e3 / (a * a * a);
         //
         // Now solve for x:
         //
         T u = exp(mu * log(nu) - nu - (1 + mu) * log(1 + mu) + mu);
         T alpha = mu;
         T u_2 = u * u;
         T u_3 = u_2 * u;
         T u_4 = u_2 * u_2;
         T u_5 = u_3 * u_2;

         T mult = u * u / 2;
         unsigned k = 2;
         T x = u + alpha * u;
         while(k < 10)
         {
            ++k;
            mult *= u / k;
            T s = mult;
            for(unsigned j = 0; j < k - 1; ++j)
            {
               s *= k * alpha + j;
            }
            x += s;
         }
         //
         // This approximation to x so far is straight out of Temme's paper
         // but isn't all that accurate, need to clean it up with some
         // Newton-Raphson iterations to stand any chance of being close
         // to the incomplete Beta inverse.  This is still *much* cheaper
         // than iterating on the full incomplete beta:
         //
         //x = 1 - x;
         u = nu - mu * log(nu) + (1 + mu) * log(1 + mu) - mu;
         x = tools::newton_raphson_iterate(
            temme_root_finder2<T>(u, alpha), x, T(0), T(1), (2 * tools::digits(x)) / 3);
         //if((nu < mu) && (x < 1 / (1 + mu)))
         //   x = 1 - x;
#ifdef BOOST_INSTRUMENT
         T x2 = tools::newton_raphson_iterate(
            temme_root_finder2<T>(u, alpha), 0.001, T(0), T(1), tools::digits(x) / 2);
         T x3 = tools::newton_raphson_iterate(
            temme_root_finder2<T>(u, alpha), 0.999, T(0), T(1), tools::digits(x) / 2);
         std::cout << "mu = " << mu << " nu = " << nu << " u = " << exp(u) << std::endl;
         std::cout << "Estimating x with Temme method 3: " << x << " (x2 = " << x2 << " x3 = " << x3 << ")" << std::endl;
#endif
         return x;
      }
   }

   T ap1 = a + 1;
   T bm1 = b - 1;
#if 0
   T lbase = boost::math::lgamma(a) + boost::math::lgamma(b) - boost::math::lgamma(a + b) + log(a) + log(z);
   T result = lbase = exp(lbase / a);
#endif
   T result = pow(a * z * boost::math::beta(a, b), 1/a);
   T lbase = result;
   result += (bm1 / ap1) * lbase * lbase;
   result += (bm1 * (a * a + 3 * b * a - a + 5 * b - 4)) / (2 * ap1 * ap1 * (a + 2)) * lbase * lbase * lbase;
#ifdef BOOST_INSTRUMENT
   std::cout << "Estimating x with Mathworld Asymptotic: " << result << std::endl;
#endif

   // clean up just in case:
   if(result > 1)
      result = 0.9999F;
   if(result < 0)
      result = 0.00000005F;

   return result;
}

template <class T, class L>
struct ibeta_roots
{
   ibeta_roots(T _a, T _b, T t, bool inv = false) 
      : a(_a), b(_b), target(t), invert(inv) {}

   std::tr1::tuple<T, T, T> operator()(const T& x)
   {
      T f = ibeta_imp(a, b, x, L(), invert, true) - target;
      T f1 = invert ? 
               -ibeta_power_terms(b, a, 1 - x, x, L(), true) 
               : ibeta_power_terms(a, b, x, 1 - x, L(), true);
      T y = 1 - x;
      if(y == 0)
         y = tools::min_value(x) * 8;
      f1 /= y * x;
      T f2 = f1 * (-y * a + (b - 2) * x + 1) / (y * x);
      if(invert)
         f2 = -f2;

      // make sure we don't have a zero derivative:
      if(f1 == 0)
         f1 = (invert ? -1 : 1) * tools::min_value(f1) * 64;

      return std::tr1::make_tuple(f, f1, f2);
   }
private:
   T a, b, target;
   bool invert;
};

template <class T, class L>
T inverse_ibeta_imp(T a, T b, T z, T guess, bool invert, int bits, L const &)
{
   //
   // Handle special cases first:
   //
   if(z == 1)
      return invert ? 0 : 1;
   if(z == 0)
      return invert ? 1 : 0;
   //
   // We need a good estimate of the error in the incomplete beta function
   // so that we don't set the desired precision too high.  Assume that 3-bits
   // are lost each time the arguments increase by a factor of 10:
   //
   using namespace std;
   int bits_lost = static_cast<int>(tools::real_cast<float>(ceil(log10((std::max)(a, b))) * 3));
   if(bits_lost < 0)
      bits_lost = 3;
   else
      bits_lost += 3;
   int precision = (std::min)(tools::digits(z) - bits_lost, bits);

   T min = 0;
   T max = 1;
   T x;
   if(z > 0.5)
      x = tools::halley_iterate(ibeta_roots<T, L>(a, b, 1 - z, !invert), guess, min, max, precision);
   else
      x = tools::halley_iterate(ibeta_roots<T, L>(a, b, z, invert), guess, min, max, precision);
#ifdef BOOST_INSTRUMENT
   std::cout << "Found x: " << x << std::endl;
   std::cout << "Forward Error: " << fabs(((invert ? ibetac(a, b, guess) : ibeta(a, b, guess)) - z) / z) << std::endl;
   if(fabs(x - guess) > 0.1)
      std::cout << "Guess was inaccurate... " << std::endl;
#endif
   return x;
}

} // namespace detail

template <class T>
T ibeta_inv(T a, T b, T x)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   T guess;
   if(x > 0.5)
      guess = 1 - detail::estimate_inverse_beta(b, a, 1 - x);
   else
      guess = detail::estimate_inverse_beta(a, b, x);
   if((guess == 0) && (x != 0))
      guess = 0.00005f;
   return static_cast<T>(detail::inverse_ibeta_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(x), static_cast<value_type>(guess), false, tools::digits(x), evaluation_type()));
}

template <class T>
T ibetac_inv(T a, T b, T x)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;
   T guess;
   if(x > 0.5)
      guess = detail::estimate_inverse_beta(a, b, 1 - x);
   else
      guess = 1 - detail::estimate_inverse_beta(b, a, x);
   if((guess == 1) && (x != 0))
      guess = 0.99F;
   return static_cast<T>(detail::inverse_ibeta_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(x), static_cast<value_type>(guess), true, tools::digits(x), evaluation_type()));
}

}} // namespaces

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_IGAMMA_INVERSE_HPP


