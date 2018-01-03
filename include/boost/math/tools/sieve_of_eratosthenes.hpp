/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 */
#ifndef BOOST_MATH_TOOLS_SIEVE_OF_ERATOSTHENESE_HPP
#define BOOST_MATH_TOOLS_SIEVE_OF_ERATOSTHENESE_HPP
#include <boost/math/tools/floor_sqrt.hpp>

namespace boost { namespace math {

template<class Z>
class sieve_of_eratosthenes {
public:
      sieve_of_eratosthenes(Z max_j) : m_max_j{max_j}
      {
          std::vector<bool> is_prime((size_t) max_j+1, true);
          is_prime[0] = false;
          is_prime[1] = false;
          Z p = 2;
          Z root_j = floor_sqrt(max_j);
          while (p <= root_j)
          {
              Z i = p + p;
              while (i <= max_j)
              {
                  is_prime[(size_t) i] = false;
                  i += p;
              }
              i = p + 1;
              while (i <= root_j && is_prime[(size_t) i] ==  false)
              {
                  ++i;
              }
              p = i;
          }
          // TODO: Tighten this up using the prime counting function.
          m_prime_list.reserve((size_t)max_j);
          for(size_t i = 0; i < is_prime.size(); ++i)
          {
              if(is_prime[i])
              {
                  m_prime_list.push_back(i);
              }
          }
          m_prime_list.shrink_to_fit();
      }

      bool is_prime(Z i) const
      {
          if (std::binary_search(m_prime_list.begin(), m_prime_list.end(), i))
          {
              return true;
          }
          if (i > m_max_j)
          {
              throw std::domain_error("Cannot determine if a number bigger than the prime table is prime.\n");
          }
          return false;
      }

      size_t prime_count() const
      {
          return m_prime_list.size();
      }

      Z prime(size_t i) const
      {
          return m_prime_list[i];
      }

private:
    std::vector<Z> m_prime_list;
    Z m_max_j;
};

}}
#endif
