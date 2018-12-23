//  (C) Copyright Anton Bikineev 2014
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_RECURRENCE_HPP_
#define BOOST_MATH_TOOLS_RECURRENCE_HPP_

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/tools/fraction.hpp>


namespace boost {
   namespace math {
      namespace tools {
         namespace detail{

            //
            // Function ratios directly from recurrence relations:
            // H. SHINTAN, Note on Miller’s recurrence algorithm, J. Sci. Hiroshima Univ. Ser. A-I
            // Math., 29 (1965), pp. 121 - 133.
            // and:
            // COMPUTATIONAL ASPECTS OF THREE-TERM RECURRENCE RELATIONS
            // WALTER GAUTSCHI
            // SIAM REVIEW Vol. 9, No. 1, January, 1967
            //
            template <class Recurrence>
            struct function_ratio_from_backwards_recurrence_fraction
            {
               typedef typename boost::remove_reference<decltype(boost::math::get<0>(std::declval<Recurrence&>()(0)))>::type value_type;
               typedef std::pair<value_type, value_type> result_type;
               function_ratio_from_backwards_recurrence_fraction(const Recurrence& r) : r(r), k(0) {}

               result_type operator()()
               {
                  value_type a, b, c;
                  boost::math::tie(a, b, c) = r(k);
                  ++k;
                  // an and bn defined as per Gauchi 1.16, not the same
                  // as the usual continued fraction a' and b's.
                  value_type bn = a / c;
                  value_type an = b / c;
                  return result_type(-bn, an);
               }

               Recurrence r;
               int k;
            };

         }  // namespace detail

         //
         // Given a stable backwards recurrence relation:
         // a f_n-1 + b f_n + c f_n+1 = 0
         // returns the ratio f_n / f_n-1
         //
         // Recurrence: a funtor that returns a tuple of the factors (a,b,c).
         // factor:     Convergence criteria, should be no less than machine epsilon.
         // max_iter:   Maximum iterations to use solving the continued fraction.
         //
         template <class Recurrence, class T>
         T function_ratio_from_backwards_recurrence(const Recurrence& r, const T& factor, boost::uintmax_t& max_iter)
         {
            detail::function_ratio_from_backwards_recurrence_fraction<Recurrence> f(r);
            return boost::math::tools::continued_fraction_a(f, factor, max_iter);
         }

         // solves usual recurrence relation for homogeneous
         // difference equation in stable forward direction
         // a(n)w(n-1) + b(n)w(n) + c(n)w(n+1) = 0
         //
         // Params:
         // get_coefs: functor returning a tuple, where
         //            get<0>() is a(n); get<1>() is b(n); get<2>() is c(n);
         // last_index: index N to be found;
         // first: w(-1);
         // second: w(0);
         //
         template <class T, class NextCoefs>
         inline T apply_recurrence_relation_forward(NextCoefs& get_coefs, unsigned last_index, T first, T second)
         {
            using std::swap;
            using boost::math::tuple;
            using boost::math::get;

            T third = 0;
            T a, b, c;

            for (unsigned k = 0; k < last_index; ++k)
            {
               tie(a, b, c) = get_coefs(k);
               // scale each part seperately to avoid spurious overflow:
               third = (a / -c) * first + (b / -c) * second;

               swap(first, second);
               swap(second, third);
            }
            
            return second;
         }

         // solves usual recurrence relation for homogeneous
         // difference equation in stable backward direction
         // a(n)w(n-1) + b(n)w(n) + c(n)w(n+1) = 0
         //
         // Params:
         // get_coefs: functor returning a tuple, where
         //            get<0>() is a(n); get<1>() is b(n); get<2>() is c(n);
         // last_index: index N to be found;
         // first: w(1);
         // second: w(0);
         //
         template <class T, class NextCoefs>
         inline T apply_recurrence_relation_backward(NextCoefs& get_coefs, unsigned last_index, T first, T second)
         {
            using std::swap;
            using boost::math::tuple;
            using boost::math::get;

            T next = 0;
            T a, b, c;

            for (unsigned k = 0; k < last_index; ++k)
            {
               tie(a, b, c) = get_coefs(-static_cast<int>(k));
               // scale each part seperately to avoid spurious overflow:
               next = (b / -a) * second + (c / -a) * first;

               swap(first, second);
               swap(second, next);
            }
            
            return second;
         }

      }
   }
} // namespaces

#endif // BOOST_MATH_TOOLS_RECURRENCE_HPP_
