//  (C) Copyright Anton Bikineev 2014
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_RECURRENCE_HPP_
#define BOOST_MATH_TOOLS_RECURRENCE_HPP_

#include <vector>
#include <algorithm>
#include <functional>

#include <boost/mpl/bool.hpp>
#include <boost/bind.hpp>

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/tools/tuple.hpp>
#include <boost/math/special_functions/next.hpp>

// There are many ways to optimize checking for close overflow
// in Olver's method. log(abs(...))-method doesn't seem to work
// fast for non-builtin types. Also it is possible to make
// file independent on boost::bind library. So we will leave some
// optimisations for future.

namespace boost {
   namespace math {
      namespace tools {

         namespace detail {

            template <class T>
            struct is_homogeneous : boost::mpl::bool_<
               boost::math::tuple_size<typename T::result_type>::value == 3u>::type {};

            // solves recurrence relation defined by homogeneous difference equation
            // a(n)w(n+1) - b(n)w(n) + c(n)w(n-1) = 0
            template <class Functor, class U, class T>
            inline std::pair<T, unsigned> olver_checked_recurrence_imp(Functor& get_coefs, const U& factor, const T& init_value, unsigned init_pos, unsigned index, boost::mpl::true_)
            {
               BOOST_MATH_STD_USING // fabs
                  using boost::math::get;
               typedef typename Functor::result_type coef_tuple;

               std::vector<T> p, e, check_ns;
               p.reserve(index); e.reserve(index); check_ns.reserve(index);

               // initialization
               coef_tuple coefs = get_coefs(init_pos + 1);
               T an = get<0>(coefs),
                  bn = get<1>(coefs),
                  cn = get<2>(coefs);

               p.push_back(0); p.push_back(1);
               e.push_back(init_value); e.push_back((cn * init_value) / an);

               T check_n = 0;
               T min_check_n = boost::math::tools::max_value<T>();
               unsigned i = 2;

               // forward recurrence
               do {
                  const T next_p = ((bn * p[i - 1]) - (cn * p[i - 2])) / an;

                  coefs = get_coefs(init_pos + i);
                  an = get<0>(coefs);
                  bn = get<1>(coefs);
                  cn = get<2>(coefs);

                  const T next_e = (cn * e[i - 1]) / an;

                  if (log(fabs((bn * p[i - 1]) - (cn * p[i - 2]))) >= log(boost::math::tools::max_value<T>()) + log(fabs(an)) ||
                     log(fabs(cn)) + log(fabs(e[i - 1])) >= log(boost::math::tools::max_value<T>()) + log(fabs(an))) // TODO: this check takes quite long time
                  {
                     typename std::vector<T>::iterator min_check_it =
                        std::find_if(check_ns.begin(), check_ns.end(), boost::bind(std::less<T>(), _1, check_ns.back() / factor)); // TODO: find better algorithm with no boost::bind
                     index = std::distance(check_ns.begin(), min_check_it);

                     break;
                  }

                  p.push_back(next_p); e.push_back(next_e);

                  check_n = fabs(e[i - 1] / (p[i - 1] * p[i]));
                  if ((i <= index) && (check_n < min_check_n))
                     min_check_n = check_n;

                  check_ns.push_back(check_n);

                  ++i;
               } while (check_n > fabs(factor * min_check_n));

               std::vector<T> w;
               w.resize(p.size());

               unsigned k = w.size();
               w[--k] = 0;

               // backward recurrence
               for (; k >= index; --k)
                  w[k - 1] = (p[k - 1] * w[k] + e[k - 1]) / p[k];

               return std::make_pair(w[index], index);
            }

            // solves recurrence relation defined by homogeneous difference equation
            // a(n)w(n+1) - b(n)w(n) + c(n)w(n-1) = d(n)
            template <class Functor, class U, class T>
            inline std::pair<T, unsigned> olver_checked_recurrence_imp(Functor& get_coefs, const U& factor, const T& init_value, unsigned init_pos, unsigned index, boost::mpl::false_)
            {
               BOOST_MATH_STD_USING // fabs
                  using boost::math::get;
               typedef typename Functor::result_type coef_tuple;

               std::vector<T> p, e, check_ns;
               p.reserve(index); e.reserve(index); check_ns.reserve(index);

               // initialization
               coef_tuple coefs = get_coefs(init_pos + 1);
               T an = get<0>(coefs),
                  bn = get<1>(coefs),
                  cn = get<2>(coefs),
                  dn = get<3>(coefs);

               p.push_back(0); p.push_back(1);
               e.push_back(init_value); e.push_back(((cn * init_value) - dn) / an);

               T check_n = 0;
               T min_check_n = boost::math::tools::max_value<T>();
               unsigned i = 2;

               // forward recurrence
               do {
                  const T next_p = ((bn * p[i - 1]) - (cn * p[i - 2]) + dn) / an;

                  coefs = get_coefs(init_pos + i);
                  an = get<0>(coefs);
                  bn = get<1>(coefs);
                  cn = get<2>(coefs);
                  dn = get<3>(coefs);

                  const T next_e = ((cn * e[i - 1]) - (dn * next_p)) / an;

                  if (log(fabs((bn * p[i - 1]) - (cn * p[i - 2]) + dn)) >= log(boost::math::tools::max_value<T>()) + log(fabs(an)) ||
                     log(fabs((cn * e[i - 1]) - (dn * next_p))) >= log(boost::math::tools::max_value<T>()) + log(fabs(an))) // TODO: this check takes quite long time
                  {
                     typename std::vector<T>::iterator min_check_it =
                        std::find_if(check_ns.begin(), check_ns.end(), boost::bind(std::less<T>(), _1, check_ns.back() / factor)); // TODO: find better algorithm with no boost::bind
                     index = std::distance(check_ns.begin(), min_check_it);

                     break;
                  }

                  p.push_back(next_p); e.push_back(next_e);

                  check_n = fabs(e[i - 1] / (p[i - 1] * p[i]));
                  if ((i <= index) && (check_n < min_check_n))
                     min_check_n = check_n;

                  check_ns.push_back(check_n);

                  ++i;
               } while (check_n > fabs(factor * min_check_n));

               std::vector<T> w;
               w.resize(p.size());

               unsigned k = w.size();
               w[--k] = 0;

               // backward recurrence
               for (; k >= index; --k)
                  w[k - 1] = (p[k - 1] * w[k] + e[k - 1]) / p[k];

               return std::make_pair(w[index], index);
            }

            // this wrapper-implementation protects us from possible overflow
            template <class Functor, class U, class T, class IsHomogeneous>
            inline T solve_recurrence_relation_by_olver_imp(Functor& get_coefs, const U& factor, unsigned index, T init_value, const IsHomogeneous& is_homogeneous)
            {
               unsigned init_pos = 0u;
               unsigned new_index = index;

               do
               {
                  std::pair<T, unsigned> result = detail::olver_checked_recurrence_imp(get_coefs, factor, init_value, init_pos, new_index, is_homogeneous);

                  init_pos += result.second;
                  new_index -= result.second;
                  init_value = result.first;

               } while (init_pos < index);

               return init_value;
            }

         } // namespace detail

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
         inline T solve_recurrence_relation_forward(NextCoefs& get_coefs, unsigned last_index, T first, T second)
         {
            using std::swap;
            using boost::math::tuple;
            using boost::math::get;

            T third = 0;
            T a, b, c;

            for (unsigned k = 0; k < last_index; ++k)
            {
               tie(a, b, c) = get_coefs(k);

               third = (a * first + b * second) / -c;

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
         inline T solve_recurrence_relation_backward(NextCoefs& get_coefs, unsigned last_index, T first, T second)
         {
            using std::swap;
            using boost::math::tuple;
            using boost::math::get;

            T next = 0;
            T a, b, c;

            for (unsigned k = 0; k < last_index; ++k)
            {
               tie(a, b, c) = get_coefs(-static_cast<int>(k));

               next = (b * second + c * first) / -a;

               swap(first, second);
               swap(second, next);
            }
            
            return second;
         }

         // solves difference equations of the following form in unstable directions:
         // a(n)w(n+1) - b(n)w(n) + c(n)w(n-1) = d(n) - inhomogeneous case
         // a(n)w(n+1) - b(n)w(n) + c(n)w(n-1) = 0    - homogeneous case
         //
         // This implementation uses Olver's algorithm because of some
         // valuable advantages as opposed to usual Miller's algorithm
         template <class Coefficients, class U, class T>
         inline T solve_recurrence_relation_by_olver(Coefficients& coefs, const U& factor, unsigned index, const T& init_value)
         {
            typedef typename detail::is_homogeneous<Coefficients>::type is_homogeneous;

            return detail::solve_recurrence_relation_by_olver_imp(coefs, factor, index, init_value, is_homogeneous());
         }

      }
   }
} // namespaces

#endif // BOOST_MATH_TOOLS_RECURRENCE_HPP_
