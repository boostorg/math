
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 Nikhar Agrawal
//  Copyright 2013 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2013 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_POLYGAMMA_DETAIL_2013_07_30_HPP_
  #define _BOOST_POLYGAMMA_DETAIL_2013_07_30_HPP_

  #include <cmath>
  #include <limits>
  #include <boost/cstdint.hpp>
  #include <boost/math/policies/policy.hpp>
  #include <boost/math/special_functions/bernoulli.hpp>
  #include <boost/math/special_functions/trunc.hpp>
  #include <boost/math/special_functions/zeta.hpp>
  #include <boost/math/special_functions/digamma.hpp>
  #include <boost/math/special_functions/sin_pi.hpp>
  #include <boost/math/special_functions/cos_pi.hpp>
  #include <boost/math/special_functions/pow.hpp>
  #include <boost/mpl/if.hpp>
  #include <boost/mpl/int.hpp>
  #include <boost/static_assert.hpp>
  #include <boost/type_traits/is_convertible.hpp>

  namespace boost { namespace math { namespace detail{

  template<class T, class Policy>
  T polygamma_atinfinityplus(const int n, const T& x, const Policy& pol, const char* function) // for large values of x such as for x> 400
  {
     // See http://functions.wolfram.com/GammaBetaErf/PolyGamma2/06/02/0001/
     BOOST_MATH_STD_USING
     //
     // sum       == current value of accumulated sum.
     // term      == value of current term to be added to sum.
     // part_term == value of current term excluding the Bernoulli number part
     //
     T term, sum, part_term;
     T x_squared = x * x;
     //
     // Start by setting part_term to:
     //
     // (n-1)! / x^(n+1)
     //
     // which is common to both the first term of the series (with k = 1)
     // and to the leading part.  
     // We can then get to the leading term by:
     //
     // part_term * (n + 2 * x) / 2
     //
     // and to the first term in the series 
     // (excluding the Bernoulli number) by:
     //
     // part_term n * (n + 1) / (2x)
     //
     // If either the factorial would overflow,
     // or the power term underflows, this just gets set to 0 and then we
     // know that we have to use logs for the initial terms:
     //
     part_term = ((n > boost::math::max_factorial<T>::value) && (n * n > tools::log_max_value<T>())) 
        ? T(0) : static_cast<T>(boost::math::factorial<T>(n - 1, pol) * pow(x, -n - 1));
     if(part_term == 0)
     {
        // Either n is very large, or the power term underflows,
        // set the initial values of part_term, term and sum via logs:
        part_term = boost::math::lgamma(n, pol) - (n + 1) * log(x);
        sum = exp(part_term + log(n + 2 * x) - boost::math::constants::ln_two<T>());
        part_term += log(T(n) * (n + 1)) - boost::math::constants::ln_two<T>() - log(x);
        part_term = exp(part_term);
     }
     else
     {
        sum = part_term * (n + 2 * x) / 2;
        part_term *= n * (n + 1) / 2;
        part_term /= x;
     }
     //
     // If the leading term is 0, so is the result:
     //
     if(sum == 0)
        return sum;

     for(unsigned k = 1;;)
     {
        term = part_term * boost::math::bernoulli_b2n<T>(k, pol);
        sum += term;
        //
        // Normal termination condition:
        //
        if(fabs(term / sum) < tools::epsilon<T>())
           break;
        //
        // Increment our counter, and move part_term on to the next value:
        //
        ++k;
        part_term *= (n + 2 * k - 2) * (n - 1 + 2 * k);
        part_term /= (2 * k - 1) * 2 * k;
        part_term /= x_squared;
        //
        // Emergency get out termination condition:
        //
        if(k > policies::get_max_series_iterations<Policy>())
        {
           return policies::raise_evaluation_error(function, "Series did not converge, closest value was %1%", sum, pol);
           break;
        }
     }
     
     if((n - 1) & 1)
        sum = -sum;

     return sum;
  }

  template<class T, class Policy>
  T polygamma_attransitionplus(const int n, const T& x, const Policy& pol, const char* function)
  {
    // See: http://functions.wolfram.com/GammaBetaErf/PolyGamma2/16/01/01/0017/

    // Use N = (0.4 * digits) + (4 * n) for target value for x:
    BOOST_MATH_STD_USING
    const int d4d  = static_cast<int>(0.4F * policies::digits_base10<T, Policy>());
    const int N = d4d + (4 * n);
    const int m    = n;
    const int iter = N - itrunc(x);

    const int minus_m_minus_one = -m - 1;

    T z(x);
    T sum0(0);
    T z_plus_k_pow_minus_m_minus_one(0);

    // Forward recursion to larger x, need to check for overflow first though:
    if(log(z + iter) * minus_m_minus_one > -tools::log_max_value<T>())
    {
       for(int k = 1; k <= iter; ++k)
       {
          z_plus_k_pow_minus_m_minus_one = pow(z, minus_m_minus_one);
          sum0 += z_plus_k_pow_minus_m_minus_one;
          z += 1;
       }
       sum0 *= boost::math::factorial<T>(n);
    }
    {
       for(int k = 1; k <= iter; ++k)
       {
          T log_term = log(z) * minus_m_minus_one + boost::math::lgamma<T>(n + 1, pol);
          sum0 += exp(log_term);
          z += 1;
       }
    }
    if((n - 1) & 1)
       sum0 = -sum0;

    return sum0 + polygamma_atinfinityplus(n, z, pol, function);
  }

  template <class T, class Policy>
  T polygamma_nearzero(int n, T x, const Policy& pol, const char* function)
  {
     BOOST_MATH_STD_USING
     //
     // If we take this expansion for polygamma: http://functions.wolfram.com/06.15.06.0003.02
     // and substitute in this expression for polygamma(n, 1): http://functions.wolfram.com/06.15.03.0009.01
     // we get an alternating series for polygamma when x is small in terms of zeta functions of
     // integer arguments (which are easy to evaluate, at least when the integer is even).
     //
     // In order to avoid spurious overflow, save the n! term for later, and rescale at the end:
     //
     T scale = boost::math::factorial<T>(n, pol);
     //
     // "factorial_part" contains everything except the zeta function
     // evaluations in each term:
     //
     T factorial_part = 1;
     //
     // "prefix" is what we'll be adding the accumulated sum to, it will
     // be n! / z^(n+1), but since we're scaling by n! it's just 
     // 1 / z^(n+1) for now:
     //
     T prefix = pow(x, n + 1);
     if(prefix == 0)
        return boost::math::policies::raise_overflow_error<T>(function, 0, pol);
     prefix = 1 / prefix;
     //
     // First term in the series is necessarily < zeta(2) < 2, so
     // ignore the sum if it will have no effect on the result anyway:
     //
     if(prefix > 2 / policies::get_epsilon<T, Policy>())
        return ((n & 1) ? 1 : -1) * 
         (tools::max_value<T>() / prefix < scale ? policies::raise_overflow_error<T>(function, 0, pol) : prefix * scale);
     //
     // Since this is an alternating sum we can accelerate convergence using
     // Algorithm 1 from "Convergence Acceleration of Alternating Series",
     // Henri Cohen, Fernando Rodriguez Villegas, and Don Zagier, 
     // Experimental Mathematics, 1999.
     // While in principle we can determine up front how many terms we will need,
     // note that often the prefix term is so large that we need no terms at all,
     // or else the series is divergent and we need rather more terms than expected.
     // The latter case we try to weed out before this function gets called, but 
     // just in case set the number of terms to an arbitrary high value, and then
     // use the usual heurists to determine when to stop, these next variables
     // correspond directly to those in Algorithm 1 of the above paper:
     //
     int nd = (int)std::min((boost::intmax_t)(boost::math::policies::digits_base10<T, Policy>() * 10), (boost::intmax_t)boost::math::policies::get_max_series_iterations<Policy>());
     T d = pow(3 + sqrt(T(8)), nd);
     d = (d + 1 / d) / 2;
     T b = -1;
     T c = -d;
     T sum = 0;
     for(int k = 0; k < nd;)
     {
        // Get the k'th term:
        T term = factorial_part * boost::math::zeta(T(k + n + 1), pol);
        // Series acceleration:
        c = b - c;
        sum = sum + c * term;
        b = (k + nd) * (k - nd) * b / ((k + 0.5) * (k + 1));
        // Termination condition:
        if(fabs(c * term) < (sum + prefix * d) * boost::math::policies::get_epsilon<T, Policy>())
           break;
        //
        // Move on k and factorial_part:
        //
        ++k;
        factorial_part *= x * (n + k) / k;
     }
     //
     // We need to add the sum to the prefix term and then
     // multiply by the scale, at each stage checking for oveflow:
     //
     sum /= d;
     if(boost::math::tools::max_value<T>() - sum < prefix)
        return boost::math::policies::raise_overflow_error<T>(function, 0, pol);
     sum += prefix;
     if(boost::math::tools::max_value<T>() / scale < sum)
        return boost::math::policies::raise_overflow_error<T>(function, 0, pol);
     sum *= scale;
     return n & 1 ? sum : -sum;
  }

  //
  // Helper function which figures out which slot our coefficient is in
  // given an angle multiplier for the cosine term of power:
  //
  template <class Table>
  typename Table::value_type::reference dereference_table(Table& table, unsigned row, unsigned power)
  {
     return table[row][power / 2];
  }



  template <class T, class Policy>
  T poly_cot_pi(int n, T x, T xc, const Policy& pol, const char* function)
  {
     // Return n'th derivative of cot(pi*x) at x, these are simply
     // tabulated for up to n = 9, beyond that it is possible to
     // calculate coefficients as follows:
     //
     // The general form of each derivative is:
     //
     // pi^n * SUM{k=0, n} C[k,n] * cos(k * pi * x) * csc^(n+1)(pi * x)
     //
     // With constant C[0,1] = -1 and all other C[k,n] = 0;
     // Then for each k < n+1:
     // C[|1 - k|, n+1]  += -(k + n + 2) * C[k, n] / 2;
     // C[k + 1, n+1]    += -(n + 2 - k) * C[k, n] / 2;
     //
     // It's worth noting however, that as well as requiring quite a bit
     // of storage space, this method has no better accuracy than recursion
     // on x to x > 0 when computing polygamma :-(
     //
     T s = fabs(x) < fabs(xc) ? boost::math::sin_pi(x) : boost::math::sin_pi(xc);
     switch(n)
     {
     case 1:
        return -constants::pi<T, Policy>() / (s * s);
     case 2:
     {
        T c = boost::math::cos_pi(x);
        return 2 * constants::pi<T, Policy>() * constants::pi<T, Policy>() * c / boost::math::pow<3>(s);
     }
     case 3:
     {
        T c = boost::math::cos_pi(2 * x);
        return -2 * boost::math::pow<3>(constants::pi<T, Policy>()) * (c + 2) / boost::math::pow<4>(s);
     }
     case 4:
     {
        T c = boost::math::cos_pi(x);
        T c2 = boost::math::cos_pi(2 * x);
        return 4 * boost::math::pow<4>(constants::pi<T, Policy>()) * (c2 + 5) * c / boost::math::pow<5>(s);
     }
     case 5:
     {
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        return -2 * boost::math::pow<5>(constants::pi<T, Policy>()) *(26 * c2 + c4 + 33) / boost::math::pow<6>(s);
     }
     case 6:
     {
        T c = boost::math::cos_pi(x);
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        return 4 * boost::math::pow<6>(constants::pi<T, Policy>()) * (56 * c2 + c4 + 123) * c / boost::math::pow<7>(s);
     }
     case 7:
     {
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        T c6 = boost::math::cos_pi(6 * x);
        return -2 * boost::math::pow<7>(constants::pi<T, Policy>()) * (1191 * c2 + 120 * c4 + c6 + 1208) / boost::math::pow<8>(s);
     }
     case 8:
     {
        T c = boost::math::cos_pi(x);
        T c3 = boost::math::cos_pi(3 * x);
        T c5 = boost::math::cos_pi(5 * x);
        T c7 = boost::math::cos_pi(7 * x);
        return 2 * boost::math::pow<8>(constants::pi<T, Policy>()) * (15619 * c + 4293 * c3 + 247 * c5 + c7) / boost::math::pow<9>(s);
     }
     case 9:
     {
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        T c6 = boost::math::cos_pi(6 * x);
        T c8 = boost::math::cos_pi(8 * x);
        return -2 * boost::math::pow<9>(constants::pi<T, Policy>()) * (88234 * c2 + 14608 * c4 + 502 * c6 + c8 + 78095) / boost::math::pow<10>(s);
     }
     case 10:
     {
        T c = boost::math::cos_pi(x);
        T c3 = boost::math::cos_pi(3 * x);
        T c5 = boost::math::cos_pi(5 * x);
        T c7 = boost::math::cos_pi(7 * x);
        T c9 = boost::math::cos_pi(9 * x);
        return 2 * boost::math::pow<10>(constants::pi<T, Policy>()) * (1310354 * c + 455192 * c3 + 47840 * c5 + 1013 * c7 + c9) / boost::math::pow<11>(s);
     }
     case 11:
     {
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        T c6 = boost::math::cos_pi(6 * x);
        T c8 = boost::math::cos_pi(8 * x);
        T c10 = boost::math::cos_pi(10 * x);
        return -2 * boost::math::pow<11>(constants::pi<T, Policy>()) * (7862124 + 9738114 * c2 + 2203488 * c4 + 152637 * c6 + 2036 * c8 + c10) / boost::math::pow<12>(s);
     }
     case 12:
     {
        T c = boost::math::cos_pi(x);
        T c3 = boost::math::cos_pi(3 * x);
        T c5 = boost::math::cos_pi(5 * x);
        T c7 = boost::math::cos_pi(7 * x);
        T c9 = boost::math::cos_pi(9 * x);
        T c11 = boost::math::cos_pi(11 * x);
        return 2 * boost::math::pow<12>(constants::pi<T, Policy>()) * (162512286 * c + 66318474 * c3 + 10187685 * c5 + 478271 * c7 + 4083 * c9 + c11) / boost::math::pow<13>(s);
     }
     case 13:
     {
        T c2 = boost::math::cos_pi(2 * x);
        T c4 = boost::math::cos_pi(4 * x);
        T c6 = boost::math::cos_pi(6 * x);
        T c8 = boost::math::cos_pi(8 * x);
        T c10 = boost::math::cos_pi(10 * x);
        T c12 = boost::math::cos_pi(12 * x);
        return -2 * boost::math::pow<13>(constants::pi<T, Policy>()) * (1137586002 + 1505621508 * c2 + 423281535 * c4 + 45533450 * c6 + 1479726 * c8 + 8178 * c10 + c12) / boost::math::pow<14>(s);
     }
     case 14:
     {
        T c = boost::math::cos_pi(x);
        T c3 = boost::math::cos_pi(3 * x);
        T c5 = boost::math::cos_pi(5 * x);
        T c7 = boost::math::cos_pi(7 * x);
        T c9 = boost::math::cos_pi(9 * x);
        T c11 = boost::math::cos_pi(11 * x);
        T c13 = boost::math::cos_pi(13 * x);
        return 2 * boost::math::pow<14>(constants::pi<T, Policy>()) * (27971176092 * c + 12843262863 * c3 + 2571742175 * c5 + 198410786 * c7 + 4537314 * c9 + 16369 * c11 + c13) / boost::math::pow<15>(s);
     }
     case 15:
     {
        return -2 * boost::math::pow<15>(constants::pi<T, Policy>()) * 
           (223769408736 + 311387598411 * boost::math::cos_pi(2 * x) 
           + 102776998928 * boost::math::cos_pi(4 * x) 
           + 15041229521 * boost::math::cos_pi(6 * x) 
           + 848090912 * boost::math::cos_pi(8 * x) 
           + 13824739 * boost::math::cos_pi(10 * x) 
           + 32752 * boost::math::cos_pi(12 * x) 
           + boost::math::cos_pi(14 * x)) / boost::math::pow<16>(s);
     }
     case 16:
     {
        return 2 * boost::math::pow<16>(constants::pi<T, Policy>()) * 
           (6382798925475 * boost::math::cos_pi(x) 
           + 3207483178157 * boost::math::cos_pi(3 * x) 
           + 782115518299 * boost::math::cos_pi(5 * x) 
           + 85383238549 * boost::math::cos_pi(7 * x) 
           + 3572085255 * boost::math::cos_pi(9 * x) 
           + 41932745 * boost::math::cos_pi(11 * x) 
           + 65519 * boost::math::cos_pi(13 * x) 
           + boost::math::cos_pi(15 * x)) / boost::math::pow<17>(s);
     }
     case 17:
     {
        return -2 * boost::math::pow<17>(constants::pi<T, Policy>()) *
           (57445190329275 + 83137223185370 * boost::math::cos_pi(2 * x)
           + 31055652948388 * boost::math::cos_pi(4 * x)
           + 5717291972382 * boost::math::cos_pi(6 * x)
           + 473353301060 * boost::math::cos_pi(8 * x)
           + 14875399450 * boost::math::cos_pi(10 * x)
           + 126781020 * boost::math::cos_pi(12 * x)
           + 131054 * boost::math::cos_pi(14 * x)
           + boost::math::cos_pi(16 * x)) / boost::math::pow<18>(s);
     }
     case 18:
     {
        return 2 * boost::math::pow<18>(constants::pi<T, Policy>()) * 
           (1865385657780650 * boost::math::cos_pi(x)
           + 1006709967915228 * boost::math::cos_pi(3 * x)
           + 285997074307300 * boost::math::cos_pi(5 * x)
           + 40457344748072 * boost::math::cos_pi(7 * x)
           + 2575022097600 * boost::math::cos_pi(9 * x)
           + 61403313100 * boost::math::cos_pi(11 * x)
           + 382439924 * boost::math::cos_pi(13 * x)
           + 262125 * boost::math::cos_pi(15 * x)
           + boost::math::cos_pi(17 * x)) / boost::math::pow<19>(s);
     }
     case 19:
     {
        return -2 * boost::math::pow<19>(constants::pi<T, Policy>()) *
           (18653856577806500 + 27862280567093358 * boost::math::cos_pi(2 * x)
           + 11485644635009424 * boost::math::cos_pi(4 * x)
           + 2527925001876036 * boost::math::cos_pi(6 * x)
           + 278794377854832 * boost::math::cos_pi(8 * x)
           + 13796160184500 * boost::math::cos_pi(10 * x)
           + 251732291184 * boost::math::cos_pi(12 * x)
           + 1151775897 * boost::math::cos_pi(14 * x)
           + 524268 * boost::math::cos_pi(16 * x)
           + boost::math::cos_pi(18 * x)) / boost::math::pow<20>(s);
     }
     case 20:
     {
        return 2 * boost::math::pow<20>(constants::pi<T, Policy>()) * 
           (679562217794156938 * boost::math::cos_pi(x)
           + 388588260723953310 * boost::math::cos_pi(3 * x)
           + 124748182104463860 * boost::math::cos_pi(5 * x)
           + 21598596303099900 * boost::math::cos_pi(7 * x)
           + 1879708669896492 * boost::math::cos_pi(9 * x)
           + 73008517581444 * boost::math::cos_pi(11 * x)
           + 1026509354985 * boost::math::cos_pi(13 * x)
           + 3464764515 * boost::math::cos_pi(15 * x)
           + 1048555 * boost::math::cos_pi(17 * x)
           + boost::math::cos_pi(19 * x)) / boost::math::pow<21>(s);
     }
     case 21:
     {
        return -2 * boost::math::pow<21>(constants::pi<T, Policy>()) *
           (7475184395735726318 + 11458681306629009100 * boost::math::cos_pi(2 * x)
           + 5119020713873609970 * boost::math::cos_pi(4 * x)
           + 1300365805079109480 * boost::math::cos_pi(6 * x)
           + 179385804170146680 * boost::math::cos_pi(8 * x)
           + 12446388300682056 * boost::math::cos_pi(10 * x)
           + 382493246941965 * boost::math::cos_pi(12 * x)
           + 4168403181210 * boost::math::cos_pi(14 * x)
           + 10414216090 * boost::math::cos_pi(16 * x)
           + 2097130 * boost::math::cos_pi(18 * x)
           + boost::math::cos_pi(20 * x)) / boost::math::pow<22>(s);
     }
     }

     //
     // We'll have to compute the corefficients up to n:
     //
#ifdef BOOST_HAS_THREADS
     static boost::detail::lightweight_mutex m;
     boost::detail::lightweight_mutex::scoped_lock l(m);
#endif
     static std::vector<std::vector<T> > table(1, std::vector<T>(1, T(-1)));

     int index = n - 1;

     if(index >= (int)table.size())
     {
        //
        // We need to compute new coefficients for the cosine terms to the derivative.
        // The following code follows immediately from differentiating
        //
        // C * cos(power * x) * csc^n(power * x)
        //
        for(int i = table.size(); i <= index; ++i)
        {
           table.push_back(std::vector<T>((i + 2) / 2, T(0)));
           for(int power = (i & 1 ? 0 : 1); power < i; power += 2)
           {
              dereference_table(table, i, std::abs(1 - power)) += -(power + i + 1) * dereference_table(table, i - 1, power) / 2;
              dereference_table(table, i, power + 1) += -(i + 1 - power) * dereference_table(table, i - 1, power) / 2;
           }
           //
           // The coefficients as calculated above grow so large so fast that we scale them all
           // by n!  And since current order = i+1 just divide each row by that as we create it:
           //
           for(unsigned j = 0; j < table[i].size(); ++j)
              table[i][j] /= (i + 1);
        }
     }

     T sum = 0;
     int power = n & 1 ? 0 : 1;
     //
     // Compute the sum of the cosine terms:
     //
     for(unsigned j = 0; j < table[index].size(); ++j)
     {
        sum += table[index][j] * boost::math::cos_pi(x * power);
        power += 2;
     }
     if(sum == 0)
        return sum;
     //
     // The remaining terms are computed using logs since the powers and factorials
     // get real large real quick:
     //
     T power_terms = n * log(boost::math::constants::pi<T>());
     if(s == 0)
        return sum * boost::math::policies::raise_overflow_error<T>(function, 0, pol);
     power_terms -= log(fabs(s)) * (n + 1);
     power_terms += boost::math::lgamma(T(n + 1));
     power_terms += log(fabs(sum));

     if(power_terms > boost::math::tools::log_max_value<T>())
        return sum * boost::math::policies::raise_overflow_error<T>(function, 0, pol);

     return exp(power_terms) * ((s < 0) && ((n + 1) & 1) ? -1 : 1) * boost::math::sign(sum);
  }

  template <class T, class Policy>
  struct polygamma_initializer
  {
     struct init
     {
        init()
        {
           // Forces initialization of our table of coefficients and mutex:
           boost::math::polygamma(30, T(-2.5), Policy());
        }
        void force_instantiate()const{}
     };
     static const init initializer;
     static void force_instantiate()
     {
        initializer.force_instantiate();
     }
  };

  template <class T, class Policy>
  const typename polygamma_initializer<T, Policy>::init polygamma_initializer<T, Policy>::initializer;
  
  template<class T, class Policy>
  inline T polygamma_imp(const int n, T x, const Policy &pol)
  {
    BOOST_MATH_STD_USING
    static const char* function = "boost::math::polygamma<%1%>(int, %1%)";
    polygamma_initializer<T, Policy>::initializer.force_instantiate();
    if(n == 0)
       return boost::math::digamma(x);
    if(n < 0)
       return policies::raise_domain_error<T>(function, "Order must be >= 0, but got %1%", n, pol);
    if(x < 0)
    {
       if(floor(x) == x)
       {
          //
          // Result is infinity if x is odd, and a pole error if x is even.
          //
          if(lltrunc(x) & 1)
             return policies::raise_overflow_error<T>(function, 0, pol);
          else
             return policies::raise_pole_error<T>(function, "Evaluation at negative integer %1%", x, pol);
       }
       //if(n < 22)
       //{
          //
          // We have tabulated the derivatives of cot(x) up to the 9th derivative, which
          // allows us to use: http://functions.wolfram.com/06.15.16.0001.01
          T z = 1 - x;
          T result = polygamma_imp(n, z, pol) + constants::pi<T, Policy>() * poly_cot_pi(n, z, x, pol, function);
          return n & 1 ? T(-result) : result;
       //}
#if 0
       //
       // Try http://functions.wolfram.com/06.15.16.0007.01
       //
       if(x <= -static_cast<int>(policies::get_max_series_iterations<Policy>()))
          return policies::raise_evaluation_error<T>(function, "Argument is outside the bounds for which we can reasonably calculate polygamma (got x = %1%)", x, pol);
       int m = boost::math::itrunc(ceil(-x));
       T z = x + m;
       T sum = 0;
       for(int k = 1; k <= m; ++k)
       {
          sum += pow(z - k, -n - 1);
       }
       sum *= boost::math::factorial<T>(n, pol);
       if(n & 1)
          sum = -sum;
       return polygamma_imp(n, z, pol) - sum;
#endif
    }
    //
    // Limit for use of small-x-series is chosen
    // so that the series doesn't go too divergent
    // in the first few terms.  Ordinarily this
    // would mean setting the limit to ~ 1 / n,
    // but since the series is accelerated, we can
    // tolerate a small amount of divergence:
    //
    T small_x_limit = std::min(T(T(5) / n), T(0.25f));
    if(x < small_x_limit)
    {
      return polygamma_nearzero(n, x, pol, function);
    }
    else if(x > 0.4F * policies::digits_base10<T, Policy>() + 4 * n)
    {
      return polygamma_atinfinityplus(n, x, pol, function);
    }
    else
    {
      return polygamma_attransitionplus(n, x, pol, function);
    }
  }

} } } // namespace boost::math::detail

#endif // _BOOST_POLYGAMMA_DETAIL_2013_07_30_HPP_

