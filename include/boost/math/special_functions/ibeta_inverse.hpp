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

//
// Helper object used by root finding
// code to convert eta to x.
//
template <class T>
struct temme_root_finder
{
   temme_root_finder(const T t_, const T a_) : t(t_), a(a_) {}

   std::tr1::tuple<T, T> operator()(const T& x)
   {
      T f = log(x) + a * log(1 - x) + t;
      T f1 = (1 / x) - (a / (1 - x));
      return std::tr1::make_tuple(f, f1);
   }
private:
   T t, a;
};
//
// See: 
// "Asymptotic Inversion of the Incomplete Beta Function"
// N.M. Temme
// Journal of Computation and Applied Mathematics 41 (1992) 145-157.
// Section 2.
//
template <class T>
T temme_method_1_ibeta_inverse(T a, T b, T z)
{
   const T r2 = sqrt(T(2));
   //
   // get the first approximation for eta from the inverse
   // error function (Eq: 2.9 and 2.10).
   //
   T eta0 = boost::math::erfc_inv(2 * z);
   eta0 /= -sqrt(a / 2);

   T terms[4] = { eta0 };
   T workspace[7];
   //
   // calculate powers:
   //
   T B = b - a;
   T B_2 = B * B;
   T B_3 = B_2 * B;
   //
   // Calculate correction terms:
   //

   // See eq following 2.15:
   workspace[0] = -B * r2 / 2;
   workspace[1] = (1 - 2 * B) / 8;
   workspace[2] = -(B * r2 / 48);
   workspace[3] = T(-1) / 192;
   workspace[4] = -B * r2 / 3840;
   terms[1] = tools::evaluate_polynomial(workspace, eta0, 5);
   // Eq Following 2.17: 
   workspace[0] = B * r2 * (3 * B - 2) / 12;
   workspace[1] = (20 * B_2 - 12 * B + 1) / 128;
   workspace[2] = B * r2 * (20 * B - 1) / 960;
   workspace[3] = (16 * B_2 + 30 * B - 15) / 4608;
   workspace[4] = B * r2 * (21 * B + 32) / 53760;
   workspace[5] = (-32 * B_2 + 63) / 368640;
   workspace[6] = -B * r2 * (120 * B + 17) / 25804480;
   terms[2] = tools::evaluate_polynomial(workspace, eta0, 7);
   // Eq Following 2.17:
   workspace[0] = B * r2 * (-75 * B_2 + 80 * B - 16) / 480;
   workspace[1] = (-1080 * B_3 + 868 * B_2 - 90 * B - 45) / 9216;
   workspace[2] = B * r2 * (-1190 * B_2 + 84 * B + 373) / 53760;
   workspace[3] = (-2240 * B_3 - 2508 * B_2 + 2100 * B - 165) / 368640;
   terms[3] = tools::evaluate_polynomial(workspace, eta0, 4);
   //
   // Bring them together to get a final estimate for eta:
   //
   T eta = tools::evaluate_polynomial(terms, 1/a, 4);
   //
   // now we need to convert eta to x, by solving the appropriate
   // quadratic equation:
   //
   T eta_2 = eta * eta;
   T c = -exp(-eta_2 / 2);
   T x;
   if(eta_2 == 0)
      x = 0.5;
   else
      x = (1 + eta * sqrt((1 + c) / eta_2)) / 2;

   BOOST_ASSERT(x >= 0);
   BOOST_ASSERT(x <= 1);
   BOOST_ASSERT(eta * (x - 0.5) >= 0);
#ifdef BOOST_INSTRUMENT
   std::cout << "Estimating x with Temme method 1: " << x << std::endl;
#endif
   return x;
}
//
// See: 
// "Asymptotic Inversion of the Incomplete Beta Function"
// N.M. Temme
// Journal of Computation and Applied Mathematics 41 (1992) 145-157.
// Section 3.
//
template <class T>
T temme_method_2_ibeta_inverse(T a, T b, T z, T r, T theta)
{
   //
   // Get first estimate for eta, see Eq 3.9 and 3.10, 
   // but note there is a typo in Eq 3.10:
   //
   T eta0 = boost::math::erfc_inv(2 * z);
   eta0 /= -sqrt(r / 2);

   T s = sin(theta);
   T c = cos(theta);
   //
   // Now we need to purturb eta0 to get eta, which we do by
   // evaluating the polynomial in 1/r at the bottom of page 151,
   // to do this we first need the error terms e1, e2 e3
   // which we'll fill into the array "terms".  Since these
   // terms are themselves polynomials, we'll need another
   // array "workspace" to calculate those...
   //
   T terms[4] = { eta0 };
   T workspace[6];
   //
   // some powers of sin(theta)cos(theta) that we'll need later:
   //
   T sc = s * c;
   T sc_2 = sc * sc;
   T sc_3 = sc_2 * sc;
   T sc_4 = sc_2 * sc_2;
   T sc_5 = sc_2 * sc_3;
   T sc_6 = sc_3 * sc_3;
   T sc_7 = sc_4 * sc_3;
   //
   // Calculate e1 and put it in terms[1], see the middle of page 151:
   //
   workspace[0] = (2 * s * s - 1) / (3 * s * c);
   static const int co1[] = { -1, -5, 5 };
   workspace[1] = -tools::evaluate_even_polynomial(co1, s, 3) / (36 * sc_2);
   static const int co2[] = { 1, 21, -69, 46 };
   workspace[2] = tools::evaluate_even_polynomial(co2, s, 4) / (1620 * sc_3);
   static const int co3[] = { 7, -2, 33, -62, 31 };
   workspace[3] = -tools::evaluate_even_polynomial(co3, s, 5) / (6480 * sc_4);
   static const int co4[] = { 25, -52, -17, 88, -115, 46 };
   workspace[4] = tools::evaluate_even_polynomial(co4, s, 6) / (90720 * sc_5);
   terms[1] = tools::evaluate_polynomial(workspace, eta0, 5);
   //
   // Now evaluate e2 and put it in terms[2]:
   //
   static const int co5[] = { 7, 12, -78, 52 };
   workspace[0] = -tools::evaluate_even_polynomial(co5, s, 4) / (405 * sc_3);
   static const int co6[] = { -7, 2, 183, -370, 185 };
   workspace[1] = tools::evaluate_even_polynomial(co6, s, 5) / (2592 * sc_4);
   static const int co7[] = { -533, 776, -1835, 10240, -13525, 5410 };
   workspace[2] = -tools::evaluate_even_polynomial(co7, s, 6) / (204120 * sc_5);
   static const int co8[] = { -1579, 3747, -3372, -15821, 45588, -45213, 15071 };
   workspace[3] = -tools::evaluate_even_polynomial(co7, s, 6) / (2099520 * sc_6);
   terms[2] = tools::evaluate_polynomial(workspace, eta0, 4);
   //
   // And e3, and put it in terms[3]:
   //
   static const int co9[] = {449, -1259, -769, 6686, -9260, 3704 };
   workspace[0] = tools::evaluate_even_polynomial(co9, s, 6) / (102060 * sc_5);
   static const int co10[] = { 63149, -151557, 140052, -727469, 2239932, -2251437, 750479 };
   workspace[1] = -tools::evaluate_even_polynomial(co10, s, 7) / (20995200 * sc_6);
   static const int co11[] = { 29233, -78755, 105222, 146879, -1602610, 3195183, -2554139, 729754 };
   workspace[2] = tools::evaluate_even_polynomial(co11, s, 8) / (36741600 * sc_7);
   terms[3] = tools::evaluate_polynomial(workspace, eta0, 3);
   //
   // Bring the correction terms together to evaluate eta, 
   // this is the last equation on page 151:
   //
   T eta = tools::evaluate_polynomial(terms, 1/r, 4);
   //
   // Now that we have eta we need to back solve for x,
   // we seek the value of x that gives eta in Eq 3.2.
   // The two methods used are described in section 5.
   //
   // Begin by defining a few variables we'll need later:
   //
   T x;
   T s_2 = s * s;
   T c_2 = c * c;
   T alpha = c / s;
   alpha *= alpha;
   T lu = (-(eta * eta) / (2 * s_2) + log(s_2) + c_2 * log(c_2) / s_2);
   //
   // Temme doesn't specify what value to switch on here, 
   // but this seems to work pretty well:
   //
   if(fabs(eta) < 0.7)
   {
      //
      // Small eta use the expansion Temme gives in the second equation
      // of section 5, it's a polynomial in eta:
      //
      workspace[0] = s * s;
      workspace[1] = s * c;
      workspace[2] = (1 - 2 * workspace[0]) / 3;
      static const int co3[] = { 1, -13, 13 };
      workspace[3] = tools::evaluate_polynomial(co3, workspace[0], 3) / (36 * s * c);
      static const int co4[] = { 1, 21, -69, 46 };
      workspace[4] = tools::evaluate_polynomial(co4, workspace[0], 4) / (270 * workspace[0] * c * c);
      x = tools::evaluate_polynomial(workspace, eta, 5);
#ifdef BOOST_INSTRUMENT
      std::cout << "Estimating x with Temme method 2 (small eta): " << x << std::endl;
#endif
   }
   else
   {
      //
      // If eta is large we need to solve Eq 3.2 more directly,
      // begin by getting an initial approximation for x from
      // the last equation on page 155, this is a polynomial in u:
      //
      T u = exp(lu);
      workspace[0] = u;
      workspace[1] = alpha;
      workspace[2] = 0;
      workspace[3] = 3 * alpha * (3 * alpha + 1) / 6;
      workspace[4] = 4 * alpha * (4 * alpha + 1) * (4 * alpha + 2) / 24;
      workspace[5] = 5 * alpha * (5 * alpha + 1) * (5 * alpha + 2) * (5 * alpha + 3) / 120;
      x = tools::evaluate_polynomial(workspace, u, 6);
      //
      // At this point we may or may not have the right answer, Eq-3.2 has
      // two solutions for x for any given eta, however the mapping in 3.2
      // is 1:1 with the sign of eta and x-sin^2(theta) being the same.
      // So we can check if we have the right root of 3.2, and if not
      // switch x for 1-x.  This transformation is motivated by the fact
      // that the distribution is *almost* symetric so 1-x will be in the right
      // ball park for the solution:
      //
      if((x - s_2) * eta < 0)
         x = 1 - x;
#ifdef BOOST_INSTRUMENT
      std::cout << "Estimating x with Temme method 2 (large eta): " << x << std::endl;
#endif
   }
   //
   // The final step is a few Newton-Raphson iterations to
   // clean up our approximation for x, this is pretty cheap
   // in general, and very cheap compared to an incomplete beta
   // evaluation.  The limits set on x come from the observation
   // that the sign of eta and x-sin^2(theta) are the same.
   //
   T lower, upper;
   if(eta < 0)
   {
      lower = 0;
      upper = s_2;
   }
   else
   {
      lower = s_2;
      upper = 1;
   }
   //
   // If our initial approximation is out of bounds then bisect:
   //
   if((x < lower) || (x > upper))
      x = (lower+upper) / 2;
   //
   // And iterate:
   //
   x = tools::newton_raphson_iterate(
      temme_root_finder<T>(-lu, alpha), x, lower, upper, (2 * tools::digits(x)) / 3);

   return x;
}

template <class T>
T temme_method_3_ibeta_inverse(T a, T b, T z)
{
   T eta0 = boost::math::gamma_Q_inv(b, z);
   eta0 /= a;
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
   T d = eta0 - mu;
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

   T eta = eta0 + e1 / a + e2 / (a * a) + e3 / (a * a * a);
   //
   // Now solve for x:
   //
   T u = eta - mu * log(eta) + (1 + mu) * log(1 + mu) - mu;
   T cross = 1 / (1 + mu);
   T lower = eta < mu ? cross : 0;
   T upper = eta < mu ? 1 : cross;
   T x = (lower + upper) / 2;
   x = tools::newton_raphson_iterate(
      temme_root_finder<T>(u, mu), x, lower, upper, (2 * tools::digits(x)) / 3);
#ifdef BOOST_INSTRUMENT
   T x2 = tools::newton_raphson_iterate(
      temme_root_finder<T>(u, mu), 0.001, T(0), T(1), tools::digits(x) / 2);
   T x3 = tools::newton_raphson_iterate(
      temme_root_finder<T>(u, mu), 0.999, T(0), T(1), tools::digits(x) / 2);
   std::cout << "mu = " << mu << " eta = " << eta << " u = " << exp(u) << std::endl;
   std::cout << "Estimating x with Temme method 3: " << x << " (x2 = " << x2 << " x3 = " << x3 << ")" << std::endl;
#endif
   return x;
}

template <class T>
T estimate_mathworld_inverse_beta(T a, T b, T p)
{
   T ap1 = a + 1;
   T bm1 = b - 1;
   T result = pow(a * p * boost::math::beta(a, b), 1/a);
   T lbase = result;
   result += (bm1 / ap1) * lbase * lbase;
   result += (bm1 * (a * a + 3 * b * a - a + 5 * b - 4)) / (2 * ap1 * ap1 * (a + 2)) * lbase * lbase * lbase;
#ifdef BOOST_INSTRUMENT
   std::cout << "Estimating x with Mathworld Asymptotic: " << result << std::endl;
#endif
   return result;
}

template <class T>
T estimate_inverse_beta(T a, T b, T p)
{
   using namespace std;
   if(p == 1)
      return 1;
   if(p == 0)
      return 0;

   T x;

   if(a + b > 5)
   {
      T minv = (std::min)(a, b);
      T maxv = (std::max)(a, b);
      if((sqrt(minv) > (maxv - minv)) && (minv > 5))
      {
         x = temme_method_1_ibeta_inverse(a, b, p);
      }
      T r = a + b;
      T theta = asin(sqrt(a / r));
      T s2t = minv / r;
      if((s2t >= 0.2) && (s2t <= 0.8) && (minv >= 10))
      {
         x = temme_method_2_ibeta_inverse(a, b, p, r, theta);
      }
      x = temme_method_3_ibeta_inverse(a, b, p);
   }
   else
   {
      x = estimate_mathworld_inverse_beta(a, b, p);
   }

   // clean up just in case:
   if(x > 1)
      x = 0.9999F;
   if(x < 0)
      x = 0.00000005F;

   return x;
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
   {
      invert = !invert;
      z = 1 - z;
   }
   bool swapx = false;
   if(guess > 0.9)
   {
      invert = !invert;
      guess = 1 - guess;
      std::swap(a, b);
      swapx = true;
   }
   x = tools::halley_iterate(ibeta_roots<T, L>(a, b, z, invert), guess, min, max, precision);
   if(swapx)
   {
      x = 1 - x;
   }
#ifdef BOOST_INSTRUMENT
   if(swapx)
   {
      invert = !invert;
      guess = 1 - guess;
      std::swap(a, b);
   }
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

   if((a <= 0) || (b <= 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the incomplete beta function inverse.");
   if((x < 0) || (x > 1))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Argument p outside the range [0,1] in the incomplete beta function inverse.");

   T guess;
   if(x > 0.5)
      guess = 1 - detail::estimate_inverse_beta(b, a, 1 - x);
   else
      guess = detail::estimate_inverse_beta(a, b, x);
   if((guess == 0) && (x != 0))
      guess = 0.00005f;
   return tools::checked_narrowing_cast<T>(detail::inverse_ibeta_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(x), static_cast<value_type>(guess), false, tools::digits(x), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

template <class T>
T ibetac_inv(T a, T b, T x)
{
   typedef typename lanczos::lanczos_traits<T>::value_type value_type;
   typedef typename lanczos::lanczos_traits<T>::evaluation_type evaluation_type;

   if((a <= 0) || (b <= 0))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Negative argument to the incomplete beta function inverse.");
   if((x < 0) || (x > 1))
      tools::domain_error<T>(BOOST_CURRENT_FUNCTION, "Argument p outside the range [0,1] in the incomplete beta function inverse.");

   T guess;
   if(x > 0.5)
      guess = detail::estimate_inverse_beta(a, b, 1 - x);
   else
      guess = 1 - detail::estimate_inverse_beta(b, a, x);
   if((guess == 1) && (x != 0))
      guess = 0.99F;
   return tools::checked_narrowing_cast<T>(detail::inverse_ibeta_imp(static_cast<value_type>(a), static_cast<value_type>(b), static_cast<value_type>(x), static_cast<value_type>(guess), true, tools::digits(x), evaluation_type()), BOOST_CURRENT_FUNCTION);
}

}} // namespaces

#endif // BOOST_MATH_SPECIAL_FUNCTIONS_IGAMMA_INVERSE_HPP


