/*
 *  (C) Copyright Nick Thompson 2018.
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  The Halton sequence for quasi-random number generation.
 */


#ifndef BOOST_MATH_TOOLS_HALTON_SEQUENCE_HPP
#define BOOST_MATH_TOOLS_HALTON_SEQUENCE_HPP
#include <algorithm>
#include <limits>
#include <cmath>
#include <random>
#include <boost/math/tools/van_der_corput.hpp>
#include <boost/math/tools/sieve_of_eratosthenes.hpp>

namespace boost { namespace math {

template<class Real>
class halton_sequence {
public:
    halton_sequence(size_t dimension, bool shuffle_coordinates = true, bool leaped = true) : m_dimension{dimension}
    {
        // We need to build a sieve that contains 'dimension' primes.
        // Since the prime counting function pi(x) < 1.26x/log(x), then we need dimension < 1.26x/log(x).
        // This inequality can be solved via the Lambert-W function, which is in progress.
        // For now, we'll just make a slight overestimate of the required sieve size:
        size_t l = log(dimension);
        size_t sieve_size = 2*dimension*l;
        sieve_of_eratosthenes<size_t> sieve(sieve_size);
        while (sieve.prime_count() < dimension + 1)
        {
            sieve_size *= 2;
            sieve = sieve_of_eratosthenes<size_t>(sieve_size);
        }

        m_primes.resize(dimension);
        for (size_t i = 0; i < m_primes.size(); ++i)
        {
            m_primes[i] = sieve.prime(i);
        }
        std::random_device rd;
        std::mt19937 eng(rd());
        if (shuffle_coordinates)
        {
            std::shuffle(m_primes.begin(), m_primes.end(), eng);
        }
        if (leaped)
        {
            // Randomly selecting a leap from the prime bases is recommended by Kocis and Whiten.
            // The leap must be prime,
            size_t index = eng() % dimension;
            m_leap = m_primes[index];
            m_primes[index] = sieve.prime(dimension);
        }
        else
        {
            m_leap = 1;
        }
        // The first few iterations of the Halton sequence have very poor properties,
        // so we set the iteration to large-ish number.
        m_iteration = m_leap*m_primes.back();
    }

    // Raw pointers are easier to wrap into other languages.
    void operator()(Real* v, size_t dimension)
    {
        assert(dimension == m_dimension);
        size_t count = m_iteration.fetch_add(m_leap);
        for (size_t i = 0; i < dimension; ++i)
        {
            v[i] = van_der_corput<Real, size_t>(count, m_primes[i]);
        }
    }

    template<class ForwardIterator>
    void operator()(ForwardIterator begin, ForwardIterator end)
    {
        size_t count = m_iteration.fetch_add(m_leap);
        auto it = begin;
        auto prime_it = m_primes.begin();
        size_t i = 0;
        while(it != end)
        {
            *it++ = van_der_corput<Real, size_t>(count, *prime_it++);
            ++i;
        }
        // Is this assertion necessary, or should the user be required to "do the right thing"?
        assert(i == m_dimension);
    }


private:
    std::atomic<size_t> m_iteration;
    size_t m_dimension;
    std::vector<size_t> m_primes;
    size_t m_leap;
};
}}
#endif
