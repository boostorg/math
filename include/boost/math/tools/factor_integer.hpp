// Copyright 2017 Nick Thompson
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_FACTOR_INTEGER_HPP
#define BOOST_MATH_TOOLS_FACTOR_INTEGER_HPP

#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <random>
#include <thread>
#include <atomic>
#include <algorithm>
#include <boost/optional.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <boost/integer/common_factor_rt.hpp>

namespace boost::math {

template<class Integer>
boost::optional<Integer> pollard_rho(Integer n, size_t num_threads = std::thread::hardware_concurrency())
{
    // Algorithm 5.23 of "The Joy of Factoring":
    using std::atomic;
    atomic<bool> done{false};
    Integer divisor = 0;
    atomic<Integer*> divisor_ptr{&divisor};
    auto f = [&n, &done, &divisor_ptr]()->void
    {
        boost::integer::gcd_evaluator<Integer>  gcd;
        std::random_device rd;
        auto seed = rd();
        std::mt19937_64 gen(seed);
        Integer b = gen() % (n-4) + 1;
        Integer s = gen() % (n-1);
        Integer A = s;
        Integer B = s;
        Integer g = 1;
        Integer iterations = 0;
        while (g == 1)
        {
            Integer C = 1;
            for(size_t i = 0; i < 100; ++i)
            {
                ++iterations;
                A = (A*A + b) % n;
                Integer tmp =  (B*B + b) % n;
                B = (tmp*tmp + b) % n;
                if (A > B)
                {
                    tmp = A - B;
                }
                else
                {
                    tmp = B - A;
                }
                Integer tmpC = (C*tmp) % n;
                if (tmpC != 0)
                {
                    C = tmpC;
                }
                else
                {
                    Integer p1 = gcd(tmp, n);
                    if (p1 > 1 && p1 < n)
                    {
                        C = p1;
                    }
                    break;
                }
                if (done)
                {
                    return;
                }
            }
            g = gcd(C, n);
        }
        //std::cout << "Iterations to factor " << n << " before completion: " << iterations << std::endl;
        if (g == n)
        {
            // Failure: Sad!
            return;
        }
        done = true;
        *divisor_ptr = g;
        return;
    };

    if (num_threads == 0)
    {
        num_threads = 1;
    }
    std::vector<std::thread> threads(num_threads);
    for (size_t i = 0; i < threads.size(); ++i)
    {
        threads[i] = std::thread(f);
    }
    std::for_each(threads.begin(), threads.end(),
                  std::mem_fn(&std::thread::join));

    if (*divisor_ptr == 0)
    {
        return {};
    }
    return *divisor_ptr;
}

template<class Integer>
std::vector<std::pair<Integer, Integer>> trial_division(Integer n)
{
    std::vector<std::pair<Integer, Integer>> pairs;
    for (Integer m = 0; m < boost::math::max_prime; ++m)
    {
        Integer j = 0;
        Integer pm = static_cast<Integer>(prime(m));
        while(n % pm == 0)
        {
            n = n/pm;
            ++j;
            if (n == 1)
            {
                break;
            }
        }

        if (j > 0)
        {
            pairs.emplace_back(pm, j);
        }
        if (n == 1)
        {
            break;
        }
    }
    if (n != 1)
    {
        pairs.emplace_back(n, 1);
    }
    return pairs;
}

}

#endif
