// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_DETAIL_TANH_SINH_CONSTANTS_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_TANH_SINH_CONSTANTS_HPP

#include <boost/type_traits/is_constructible.hpp>
#include <boost/lexical_cast.hpp>

namespace boost {
   namespace math {
      namespace quadrature {
         namespace detail {

            template <class Real>
            inline Real abscissa_at_t(const Real& t)
            {
               using std::tanh;
               using std::sinh;
               return tanh(constants::half_pi<Real>()*sinh(t));
            }
            template <class Real>
            inline Real weight_at_t(const Real& t)
            {
               using std::cosh;
               using std::sinh;
               Real cs = cosh(half_pi<Real>() * sinh(t));  
               return constants::half_pi<Real>() * cosh(t) / (cs * cs);
            }
            template <class Real>
            inline Real abscissa_complement_at_t(const Real& t)
            {
               using std::cosh;
               using std::exp;
               using std::sinh;
               Real u2 = constants::half_pi<Real>() * sinh(t);  
               return 1 / (exp(u2) *cosh(u2));
            }
            template <class Real>
            inline Real t_from_abscissa_complement(const Real& x)
            {
               using std::log;
               using std::sqrt;
               Real l = log(sqrt((2 - x) / x));
               return log((sqrt(4 * l * l + constants::pi<Real>() * constants::pi<Real>()) + 2 * l) / constants::pi<Real>());
            };

            template<class Real, int precision>
            class tanh_sinh_detail_constants
            {
            public:
               tanh_sinh_detail_constants(std::size_t max_refinements, std::size_t initial_commit = 4) : m_max_refinements(max_refinements), m_committed_refinements(initial_commit)
               {
                  using std::tanh;
                  using std::sinh;
                  using std::asinh;
                  using std::atanh;
                  using boost::math::constants::half_pi;
                  using boost::math::constants::pi;
                  using boost::math::constants::two_div_pi;

                  m_inital_row_length = 7;
                  std::size_t first_complement = 0;
                  m_t_max = m_inital_row_length;
                  m_t_crossover = t_from_abscissa_complement(Real(0.5f));

                  m_abscissas.assign(m_max_refinements + 1, std::vector<Real>());
                  m_weights.assign(m_max_refinements + 1, std::vector<Real>());
                  m_first_complements.assign(m_max_refinements + 1, 0);
                  //
                  // First row is special:
                  //
                  Real h = m_t_max / m_inital_row_length;

                  std::vector<Real> temp(m_inital_row_length + 1, Real(0));
                  for (std::size_t i = 0; i < m_inital_row_length; ++i)
                  {
                     Real t = h * i;
                     if ((first_complement < 0) && (t < m_t_crossover))
                        first_complement = i;
                     temp[i] = t < m_t_crossover ? abscissa_at_t(t) : -abscissa_complement_at_t(t);
                  }
                  temp[m_inital_row_length] = abscissa_complement_at_t(m_t_max);
                  m_abscissas[0].swap(temp);
                  m_first_complements[0] = first_complement;
                  temp.assign(m_inital_row_length + 1, Real(0));
                  for (std::size_t i = 0; i < m_inital_row_length; ++i)
                     temp[i] = weight_at_t(Real(h * i));
                  temp[m_inital_row_length] = weight_at_t(m_t_max);
                  m_weights[0].swap(temp);

                  for (std::size_t row = 1; row <= initial_commit; ++row)
                  {
                     h /= 2;
                     first_complement = 0;

                     for (Real pos = h; pos < m_t_max; pos += 2 * h)
                     {
                        if (pos < m_t_crossover)
                           ++first_complement;
                        temp.push_back(pos < m_t_crossover ? abscissa_at_t(pos) : -abscissa_complement_at_t(pos));
                     }
                     m_abscissas[row].swap(temp);
                     m_first_complements[row] = first_complement;
                     for (Real pos = h; pos < m_t_max; pos += 2 * h)
                        temp.push_back(weight_at_t(pos));
                     m_weights[row].swap(temp);
                  }
               }

               void extend_refinements()const
               {
                  ++m_committed_refinements;
                  std::size_t row = m_committed_refinements;
                  Real h = ldexp(Real(1), -static_cast<int>(row));
                  std::size_t first_complement = 0;
                  for (Real pos = h; pos < m_t_max; pos += 2 * h)
                  {
                     if (pos < m_t_crossover)
                        ++first_complement;
                     m_abscissas[row].push_back(pos < m_t_crossover ? abscissa_at_t(pos) : -abscissa_complement_at_t(pos));
                  }
                  m_first_complements[row] = first_complement;
                  for (Real pos = h; pos < m_t_max; pos += 2 * h)
                     m_weights[row].push_back(weight_at_t(pos));
               }

               const std::vector<Real>& get_abscissa_row(std::size_t n)const
               {
                  if (m_committed_refinements < n)
                     extend_refinements();
                  BOOST_ASSERT(m_committed_refinements >= n);
                  return m_abscissas[n];
               }
               const std::vector<Real>& get_weight_row(std::size_t n)const
               {
                  if (m_committed_refinements < n)
                     extend_refinements();
                  BOOST_ASSERT(m_committed_refinements >= n);
                  return m_weights[n];
               }
               std::size_t get_first_complement_index(std::size_t n)const
               {
                  if (m_committed_refinements < n)
                     extend_refinements();
                  BOOST_ASSERT(m_committed_refinements >= n);
                  return m_first_complements[n];
               }

            protected:
               mutable std::vector<std::vector<Real>> m_abscissas;
               mutable std::vector<std::vector<Real>> m_weights;
               mutable std::vector<std::size_t>       m_first_complements;
               std::size_t                       m_max_refinements, m_inital_row_length;
               mutable std::size_t               m_committed_refinements;
               Real m_t_max, m_t_crossover;
            };


         }
      }
   }
}

#endif


