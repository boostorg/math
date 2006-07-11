//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_NEWTON_SOLVER_HPP
#define BOOST_MATH_TOOLS_NEWTON_SOLVER_HPP

#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable: 4512)
#endif
#include <boost/tr1/tuple.hpp>
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif
#include <boost/cstdint.hpp>
#include <boost/assert.hpp>
#include <boost/throw_exception.hpp>

#include <utility>
#include <cmath>
#include <stdexcept>

namespace boost{ namespace math{ namespace tools{

namespace detail{

template <class F, class T>
void handle_zero_derivative(F f,
                            T& last_f0,
                            const T& f0,
                            T& delta,
                            T& result,
                            T& guess,
                            const T& min,
                            const T& max)
{
   if(last_f0 == 0)
   {
      // this must be the first iteration, pretend that we had a
      // previous one at either min or max:
      if(result == min)
      {
         guess = max;
      }
      else
      {
         guess = min;
      }
      last_f0 = std::tr1::get<0>(f(guess));
      delta = guess - result;
   }
   if(last_f0 * f0 < 0)
   {
      // we've crossed over so move in opposite direction to last step:
      if(delta < 0)
      {
         delta = (result - min) / 2;
      }
      else
      {
         delta = (result - max) / 2;
      }
   }
   else
   {
      // move in same direction as last step:
      if(delta < 0)
      {
         delta = (result - max) / 2;
      }
      else
      {
         delta = (result - min) / 2;
      }
   }
}

} // namespace

template <class F, class T>
T bisect(F f, T min, T max, int digits, boost::uintmax_t& max_iter)
{
   T guess = (max + min) / 2;
   T factor = static_cast<T>(ldexp(1.0, 1 - digits));

   T fmin = f(min);
   T fmax = f(max);
   T fguess = f(guess);
   if(fmin == 0)
      return min;
   if(fmax == 0)
      return max;
   if(fguess == 0)
      return guess;

   if(min >= max)
   {
      std::invalid_argument e("Arguments in wrong order in boost::math::tools::bisect");
      boost::throw_exception(e);
   }
   if(fmin * fmax >= 0)
   {
      std::logic_error e("No change of sign in boost::math::tools::bisect, either there is no root to find, or there are multiple roots in the interval.");
      boost::throw_exception(e);
   }

   boost::uintmax_t count = max_iter;

   do
   {
      if(fguess == 0)
         break;
      else if(fguess * fmin < 0)
      {
         max = guess;
         fmax = fguess;
         guess = (guess + min) / 2;
         fguess = f(guess);
         if(guess == max)
            break;
      }
      else
      {
         min = guess;
         fmin = fguess;
         guess = (guess + max) / 2;
         fguess = f(guess);
         if(guess == min)
            break;
      }
   }while(--max_iter && (max - min > factor * guess));

   max_iter -= count;

   return guess;
}

template <class F, class T>
T bisect(F f, T min, T max, int digits)
{
   boost::uintmax_t m = (std::numeric_limits<boost::uintmax_t>::max)();
   return bisect(f, min, max, digits, m);
}

template <class F, class T>
T newton_raphson_iterate(F f, T guess, T min, T max, int digits, boost::uintmax_t& max_iter)
{
   using namespace std;

   T f0(0), f1, last_f0(0);
   T result = guess;

   T factor = static_cast<T>(ldexp(1.0, 1 - digits));
   T delta = 1;
   T delta1 = tools::max_value<T>();
   T delta2 = tools::max_value<T>();

   boost::uintmax_t count(max_iter);

   do{
      last_f0 = f0;
      delta2 = delta1;
      delta1 = delta;
      std::tr1::tie(f0, f1) = f(result);
      if(0 == f0)
         break;
      if(f1 == 0)
      {
         // Oops zero derivative!!!
#ifdef BOOST_INSTRUMENT
         std::cout << "Newton iteration, zero derivative found" << std::endl;
#endif
         detail::handle_zero_derivative(f, last_f0, f0, delta, result, guess, min, max);
      }
      else
      {
         delta = f0 / f1;
      }
#ifdef BOOST_INSTRUMENT
      std::cout << "Newton iteration, delta = " << delta << std::endl;
#endif
      if(fabs(delta * 2) > fabs(delta2))
      {
         // last two steps haven't converged, try bisection:
         delta = (delta > 0) ? (result - min) / 2 : (result - max) / 2;
      }
      guess = result;
      result -= delta;
      if(result <= min)
      {
         delta = 0.5F * (guess - min);
         result = guess - delta;
         if((result == min) || (result == max))
            break;
      }
      else if(result >= max)
      {
         delta = 0.5F * (guess - max);
         result = guess - delta;
         if((result == min) || (result == max))
            break;
      }
      // update brackets:
      if(delta > 0)
         max = guess;
      else
         min = guess;
   }while(--count && (fabs(result * factor) < fabs(delta)));

   max_iter -= count;

#ifdef BOOST_INSTRUMENT
   std::cout << "Newton Raphson iteration, final count = " << max_iter << std::endl;

   static boost::uintmax_t max_count = 0;
   if(max_iter > max_count)
   {
      max_count = max_iter;
      std::cout << "Maximum iterations: " << max_iter << std::endl;
   }
#endif

   return result;
}

template <class F, class T>
T newton_raphson_iterate(F f, T guess, T min, T max, int digits)
{
   boost::uintmax_t m = (std::numeric_limits<boost::uintmax_t>::max)();
   return newton_raphson_iterate(f, guess, min, max, digits, m);
}

template <class F, class T>
T halley_iterate(F f, T guess, T min, T max, int digits, boost::uintmax_t& max_iter)
{
   using namespace std;

   T f0(0), f1, f2;
   T result = guess;

   T factor = static_cast<T>(ldexp(1.0, 1 - digits));
   T delta = (std::max)(10000000 * guess, T(10000000));  // arbitarily large delta
   T last_f0 = 0;
   T delta1 = delta;
   T delta2 = delta;

#ifdef BOOST_INSTRUMENT
   std::cout << "Halley iteration, limit = " << factor << std::endl;
#endif

   boost::uintmax_t count(max_iter);

   do{
      last_f0 = f0;
      delta2 = delta1;
      delta1 = delta;
      std::tr1::tie(f0, f1, f2) = f(result);
      if(0 == f0)
         break;
      if((f1 == 0) && (f2 == 0))
      {
         // Oops zero derivative!!!
#ifdef BOOST_INSTRUMENT
         std::cout << "Halley iteration, zero derivative found" << std::endl;
#endif
         detail::handle_zero_derivative(f, last_f0, f0, delta, result, guess, min, max);
      }
      else
      {
         if(f2 != 0)
         {
            T denom = 2 * f0;
            T num = 2 * f1 - f0 * (f2 / f1);
            if((fabs(num) < 1) && (fabs(denom) >= fabs(num) * tools::max_value<T>()))
            {
               // possible overflow, use Newton step:
               delta = f0 / f1;
            }
            else
               delta = denom / num;
            if(delta * f1 / f0 < 0)
            {
               // probably cancellation error, try a Newton step instead:
               delta = f0 / f1;
            }
         }
         else
            delta = f0 / f1;
      }
#ifdef BOOST_INSTRUMENT
      std::cout << "Halley iteration, delta = " << delta << std::endl;
#endif
      T convergence = fabs(delta / delta2);
      if((convergence > 0.5) && (convergence < 2))
      {
         // last two steps haven't converged, try bisection:
         delta = (delta > 0) ? (result - min) / 2 : (result - max) / 2;
      }
      guess = result;
      result -= delta;

      // check for out of bounds step:
      if(result <= min)
      {
         delta = (guess - min) / 2;
         result = guess - delta;
         if((result == min) || (result == max))
            break;
      }
      else if(result >= max)
      {
         delta = (guess - max) / 2;
         result = guess - delta;
         if((result == min) || (result == max))
            break;
      }
      // update brackets:
      if(delta > 0)
         max = guess;
      else
         min = guess;
   }while(--count && (fabs(result * factor) < fabs(delta)));

   max_iter -= count;

#ifdef BOOST_INSTRUMENT
   std::cout << "Halley iteration, final count = " << max_iter << std::endl;

   static boost::uintmax_t max_count = 0;
   if(max_iter > max_count)
   {
      max_count = max_iter;
      std::cout << "Maximum iterations: " << max_iter << std::endl;
   }
#endif

   return result;
}

template <class F, class T>
T halley_iterate(F f, T guess, T min, T max, int digits)
{
   boost::uintmax_t m = (std::numeric_limits<boost::uintmax_t>::max)();
   return halley_iterate(f, guess, min, max, digits, m);
}

template <class F, class T>
T schroeder_iterate(F f, T guess, T min, T max, int digits, boost::uintmax_t& max_iter)
{
   using namespace std;

   T f0(0), f1, f2, last_f0(0);
   T result = guess;

   T factor = static_cast<T>(ldexp(1.0, 1 - digits));
   T delta = 0;
   T delta1 = tools::max_value<T>();
   T delta2 = tools::max_value<T>();

#ifdef BOOST_INSTRUMENT
   std::cout << "Schroeder iteration, limit = " << factor << std::endl;
#endif

   boost::uintmax_t count(max_iter);

   do{
      last_f0 = f0;
      delta2 = delta1;
      delta1 = delta;
      std::tr1::tie(f0, f1, f2) = f(result);
      if(0 == f0)
         break;
      if((f1 == 0) && (f2 == 0))
      {
         // Oops zero derivative!!!
#ifdef BOOST_INSTRUMENT
         std::cout << "Halley iteration, zero derivative found" << std::endl;
#endif
         detail::handle_zero_derivative(f, last_f0, f0, delta, result, guess, min, max);
      }
      else
      {
         T ratio = f0 / f1;
         if(ratio / result < 0.1)
         {
            delta = ratio + (f2 / (2 * f1)) * ratio * ratio;
            // check second derivative doesn't over compensate:
            if(delta * ratio < 0)
               delta = ratio;
         }
         else
            delta = ratio;  // fall back to Newton iteration.
      }
      if(fabs(delta * 2) > fabs(delta2))
      {
         // last two steps haven't converged, try bisection:
         delta = (delta > 0) ? (result - min) / 2 : (result - max) / 2;
      }
      guess = result;
      result -= delta;
#ifdef BOOST_INSTRUMENT
      std::cout << "Halley iteration, delta = " << delta << std::endl;
#endif
      if(result <= min)
      {
         delta = 0.5F * (guess - min);
         result = guess - delta;
         if((result == min) || (result == max))
            break;
      }
      else if(result >= max)
      {
         delta = 0.5F * (guess - max);
         result = guess - delta;
         if((result == min) || (result == max))
            break;
      }
      // update brackets:
      if(delta > 0)
         max = guess;
      else
         min = guess;
   }while(--count && (fabs(result * factor) < fabs(delta)));

   max_iter -= count;

#ifdef BOOST_INSTRUMENT
   std::cout << "Schroeder iteration, final count = " << max_iter << std::endl;

   static boost::uintmax_t max_count = 0;
   if(max_iter > max_count)
   {
      max_count = max_iter;
      std::cout << "Maximum iterations: " << max_iter << std::endl;
   }
#endif

   return result;
}

template <class F, class T>
T schroeder_iterate(F f, T guess, T min, T max, int digits)
{
   boost::uintmax_t m = (std::numeric_limits<boost::uintmax_t>::max)();
   return schroeder_iterate(f, guess, min, max, digits, m);
}

} // namespace tools
} // namespace math
} // namespace boost


#endif // BOOST_MATH_TOOLS_NEWTON_SOLVER_HPP

