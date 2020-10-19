// Copyright 2020 Matt Borland
//
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_JM_HPP
#define BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_JM_HPP

#include <boost/assert.hpp>
#include <vector>
#include <iterator>
#include <cmath>
#include <thread>
#include <memory>
#include <cstring>
#include <future>
#include <numeric>
#include <algorithm>
#include <execution>
#include <type_traits>
#include <boost/type_traits/has_post_increment.hpp>
#include <boost/dynamic_bitset.hpp>

#include <intrin.h>

namespace jm { namespace detail
{
// https://mathworld.wolfram.com/SieveofEratosthenes.html
template <class I>
struct simple_bitset
{
   std::unique_ptr<I[]> bits;
   std::size_t m_size;
   simple_bitset(std::size_t n) : bits(new I[n / (sizeof(I) * CHAR_BIT) + (n % (sizeof(I) * CHAR_BIT) ? 1 : 0)]), m_size(n)
   {
      std::memset(bits.get(), 0xff, n / CHAR_BIT + (n % CHAR_BIT ? 1 : 0));
   }
   static constexpr std::size_t ln2(std::size_t n)
   {
      return n <= 1 ? 0 : 1 + ln2(n >> 1);
   }
   I test(std::size_t n)const
   {
      constexpr I masks[] = { static_cast<I>(1uLL), static_cast<I>(2uLL), static_cast<I>(4uLL), static_cast<I>(8uLL), static_cast<I>(16uLL), static_cast<I>(32uLL), static_cast<I>(64uLL), static_cast<I>(128uLL), static_cast<I>(256uLL),
         static_cast<I>(1uLL << 9), static_cast<I>(1uLL << 10), static_cast<I>(1uLL << 11), static_cast<I>(1uLL << 12), static_cast<I>(1uLL << 13), static_cast<I>(1uLL << 14), static_cast<I>(1uLL << 15), static_cast<I>(1uLL << 16),
         static_cast<I>(1uLL << 17), static_cast<I>(1uLL << 18), static_cast<I>(1uLL << 19), static_cast<I>(1uLL << 20), static_cast<I>(1uLL << 21), static_cast<I>(1uLL << 22), static_cast<I>(1uLL << 23), static_cast<I>(1uLL << 24),
         static_cast<I>(1uLL << 25), static_cast<I>(1uLL << 26), static_cast<I>(1uLL << 27), static_cast<I>(1uLL << 28), static_cast<I>(1uLL << 29), static_cast<I>(1uLL << 30), static_cast<I>(1uLL << 31), static_cast<I>(1uLL << 32),
         static_cast<I>(1uLL << 33), static_cast<I>(1uLL << 34), static_cast<I>(1uLL << 35), static_cast<I>(1uLL << 36), static_cast<I>(1uLL << 37), static_cast<I>(1uLL << 38), static_cast<I>(1uLL << 39), static_cast<I>(1uLL << 40),
         static_cast<I>(1uLL << 41), static_cast<I>(1uLL << 42), static_cast<I>(1uLL << 43), static_cast<I>(1uLL << 44), static_cast<I>(1uLL << 45), static_cast<I>(1uLL << 46), static_cast<I>(1uLL << 47), static_cast<I>(1uLL << 48),
         static_cast<I>(1uLL << 49), static_cast<I>(1uLL << 50), static_cast<I>(1uLL << 51), static_cast<I>(1uLL << 52), static_cast<I>(1uLL << 53), static_cast<I>(1uLL << 54), static_cast<I>(1uLL << 55), static_cast<I>(1uLL << 56),
         static_cast<I>(1uLL << 57), static_cast<I>(1uLL << 58), static_cast<I>(1uLL << 59), static_cast<I>(1uLL << 60), static_cast<I>(1uLL << 61), static_cast<I>(1uLL << 62), static_cast<I>(1uLL << 63),
      };
      constexpr I mask = (sizeof(I) * CHAR_BIT) - 1;
      constexpr std::size_t shift = ln2(sizeof(I) * CHAR_BIT);
      //BOOST_ASSERT((n >> shift) < (m_size / (sizeof(I) * CHAR_BIT) + (m_size % (sizeof(I) * CHAR_BIT) ? 1 : 0)));
      return bits[n >> shift] & masks[n & mask];
   }
   void clear(std::size_t n)
   {
      constexpr I masks[] = { static_cast<I>(~1uLL), static_cast<I>(~2uLL), static_cast<I>(~4uLL), static_cast<I>(~8uLL), static_cast<I>(~16uLL), static_cast<I>(~32uLL), static_cast<I>(~64uLL), static_cast<I>(~128uLL), static_cast<I>(~256uLL),
         static_cast<I>(~(1uLL << 9)), static_cast<I>(~(1uLL << 10)), static_cast<I>(~(1uLL << 11)), static_cast<I>(~(1uLL << 12)), static_cast<I>(~(1uLL << 13)), static_cast<I>(~(1uLL << 14)), static_cast<I>(~(1uLL << 15)), static_cast<I>(~(1uLL << 16)),
         static_cast<I>(~(1uLL << 17)), static_cast<I>(~(1uLL << 18)), static_cast<I>(~(1uLL << 19)), static_cast<I>(~(1uLL << 20)), static_cast<I>(~(1uLL << 21)), static_cast<I>(~(1uLL << 22)), static_cast<I>(~(1uLL << 23)), static_cast<I>(~(1uLL << 24)),
         static_cast<I>(~(1uLL << 25)), static_cast<I>(~(1uLL << 26)), static_cast<I>(~(1uLL << 27)), static_cast<I>(~(1uLL << 28)), static_cast<I>(~(1uLL << 29)), static_cast<I>(~(1uLL << 30)), static_cast<I>(~(1uLL << 31)), static_cast<I>(~(1uLL << 32)),
         static_cast<I>(~(1uLL << 33)), static_cast<I>(~(1uLL << 34)), static_cast<I>(~(1uLL << 35)), static_cast<I>(~(1uLL << 36)), static_cast<I>(~(1uLL << 37)), static_cast<I>(~(1uLL << 38)), static_cast<I>(~(1uLL << 39)), static_cast<I>(~(1uLL << 40)),
         static_cast<I>(~(1uLL << 41)), static_cast<I>(~(1uLL << 42)), static_cast<I>(~(1uLL << 43)), static_cast<I>(~(1uLL << 44)), static_cast<I>(~(1uLL << 45)), static_cast<I>(~(1uLL << 46)), static_cast<I>(~(1uLL << 47)), static_cast<I>(~(1uLL << 48)),
         static_cast<I>(~(1uLL << 49)), static_cast<I>(~(1uLL << 50)), static_cast<I>(~(1uLL << 51)), static_cast<I>(~(1uLL << 52)), static_cast<I>(~(1uLL << 53)), static_cast<I>(~(1uLL << 54)), static_cast<I>(~(1uLL << 55)), static_cast<I>(~(1uLL << 56)),
         static_cast<I>(~(1uLL << 57)), static_cast<I>(~(1uLL << 58)), static_cast<I>(~(1uLL << 59)), static_cast<I>(~(1uLL << 60)), static_cast<I>(~(1uLL << 61)), static_cast<I>(~(1uLL << 62)), static_cast<I>(~(1uLL << 63)),
      };
      constexpr I mask = (sizeof(I) * CHAR_BIT) - 1;
      constexpr std::size_t shift = ln2(sizeof(I) * CHAR_BIT);
      //BOOST_ASSERT((n >> shift) < (m_size / (sizeof(I) * CHAR_BIT) + (m_size % (sizeof(I) * CHAR_BIT) ? 1 : 0)));
      bits[n >> shift] &= masks[n & mask];
   }
   std::size_t size()const { return m_size; }
   void reset(){ std::memset(bits.get(), 0xff, m_size / CHAR_BIT + (m_size % CHAR_BIT ? 1 : 0)); }
};

template <>
inline std::uint64_t simple_bitset<std::uint64_t>::test(std::size_t n)const
{
   constexpr std::uint64_t mask = 63;
   constexpr std::size_t shift = 6;
   return _bittest64(reinterpret_cast<long long*>(bits.get()) + (n >> shift), n & mask);
}
template <>
inline void simple_bitset<std::uint64_t>::clear(std::size_t n)
{
   constexpr std::uint64_t mask = 63;
   constexpr std::size_t shift = 6;
   _bittestandreset64(reinterpret_cast<long long*>(bits.get()) + (n >> shift), n & mask);
}

struct dynamic_bitset_wrapper
{
   boost::dynamic_bitset<> data;

   bool test(std::size_t n)const
   {
      return data.test(n);
   }
   void clear(std::size_t n)
   {
      data.set(n, false);
   }

   dynamic_bitset_wrapper(std::size_t n)
   {
      data.resize(n, true);
   }
   std::size_t size()const { return data.size(); }
   void reset() 
   {
      data.set(0, data.size(), true);
   }
};

typedef simple_bitset<std::uint64_t> bitmask_type;
//typedef dynamic_bitset_wrapper bitmask_type;

template <class T>
constexpr bool has_output_iterator_terminated(const T&)
{
   return false;
}

template <class OutputIterator, class Container>
struct dual_output_iterator
{
   OutputIterator out;
   Container& container;

   dual_output_iterator(OutputIterator o, Container& c) : out(o), container(c) {}

   struct proxy
   {
      dual_output_iterator& ref;
      proxy(dual_output_iterator& o) : ref(o) {}
      ~proxy() { ++ref.out; }
      dual_output_iterator& operator*() { return ref; }
   };
   dual_output_iterator& operator++() { ++out; return *this; }
   proxy operator++(int) { return proxy(*this); }
   dual_output_iterator& operator*() { return *this; }
   template <class Value>
   dual_output_iterator& operator=(const Value& val)
   {
      *out = val;
      container.push_back(val);
      return *this;
   }
};

template<class Integer, class Sieve, class OutputIterator>
inline bool linear_sieve_classical(Sieve& masks, OutputIterator out)
{
   Integer max_sqr = Integer(1) << (std::numeric_limits<Integer>::digits / 2 - 1);
   *out++ = 2u;
   for (Integer index = 1; index < masks.size(); ++index)
   {
      if (masks.test(index))
      {
         *out++ = 2 * index + 1;
         if (has_output_iterator_terminated(out))
            return false;
         if(index < max_sqr)
            for (Integer clear = index * (2 * index + 2); clear < masks.size(); clear += index * 2 + 1)
               masks.clear(clear);
      }
   }
   return true;
}

template<class Container, class Sieve, class Integer, class OutputIterator>
inline bool linear_sieve_classical_segment(Container& primes, Sieve& masks, Integer start_offset, Integer max_points, OutputIterator out, bool output_to_container)
{
   masks.reset();
   Integer limit = static_cast<Integer>(std::sqrt(static_cast<double>(start_offset + max_points)) + 1);
   // Begin by striking out odd numbered multiples of all the primes we have so far.
   // 1-based index, we only have odd numbers in the sieve so don't need to handle 2:
   for (std::size_t i = 1; i < primes.size(); ++i)
   {
      Integer prime = primes[i];
      if (prime > limit)
         break;
      Integer index = prime - start_offset % prime;
      if ((index & 1) == 0)
         index += prime;
      index >>= 1;
      for (; index < max_points / 2; index += prime)
         masks.clear(index);
   }
   //
   // Now walk through the sieve outputting primes.
   //
   for (Integer index = 0; index < max_points / 2; ++index)
   {
      if (masks.test(index))
      {
         *out++ = start_offset + 2 * index + 1;
         if (output_to_container)
            primes.push_back(start_offset + 2 * index + 1);
         if (has_output_iterator_terminated(out))
            return false;
      }
   }
   return true;
}


// 4096 is where benchmarked performance of linear_sieve begins to diverge
template<class Integer>
const Integer linear_sieve_limit = Integer(524288); // Constexpr does not work with boost::multiprecision types

template<class Container, class Integer, class OutputIterator>
inline bool linear_sieve_classical_segment_threaded(std::atomic<Integer>* current_max_processed_value, std::mutex* lock, Container* primes, Integer start_offset, Integer end_offset, Integer stride, OutputIterator out, bool output_to_container)
{
   bitmask_type masks(linear_sieve_limit<Integer> / 2);

   std::unique_lock<std::mutex> l(*lock);
   std::size_t prime_count = primes->size();
   l.unlock();

   while (start_offset < end_offset)
   {
      Integer max_points = end_offset - start_offset > linear_sieve_limit<Integer> ? linear_sieve_limit<Integer> : end_offset - start_offset;
      Integer limit = static_cast<Integer>(std::sqrt(static_cast<double>(start_offset + max_points)) + 1);
      // Begin by striking out odd numbered multiples of all the primes we have so far.
      // 1-based index, we only have odd numbers in the sieve so don't need to handle 2:
      for (std::size_t i = 1; i < prime_count; ++i)
      {
         Integer prime = (*primes)[i];
         if (prime > limit)
            break;
         Integer index = prime - start_offset % prime;
         if ((index & 1) == 0)
            index += prime;
         index >>= 1;
         for (; index < max_points / 2; index += prime)
            masks.clear(index);
      }
      //
      // Now we must wait until the previous thread has completed the segment before this one:
      //
      while (current_max_processed_value->load() != start_offset)
         std::this_thread::yield();
      //
      // Maybe process all the primes we didn't have access to in the loop above:
      //
      if ((*primes)[prime_count - 1] < limit)
      {
         l.lock();
         prime_count = primes->size();
         l.unlock();
         // Begin by striking out odd numbered multiples of all the primes we have so far.
         // 1-based index, we only have odd numbers in the sieve so don't need to handle 2:
         for (std::size_t i = 1; i < prime_count; ++i)
         {
            Integer prime = (*primes)[i];
            if (prime > limit)
               break;
            Integer index = prime - start_offset % prime;
            if ((index & 1) == 0)
               index += prime;
            index >>= 1;
            for (; index < max_points / 2; index += prime)
               masks.clear(index);
         }
      }
      //
      // Now walk through the sieve outputting primes.
      //
      l.lock();
      for (Integer index = 0; index < max_points / 2; ++index)
      {
         if (masks.test(index))
         {
            *out++ = start_offset + 2 * index + 1;
            if (output_to_container)
               primes->push_back(start_offset + 2 * index + 1);
            if (has_output_iterator_terminated(out))
               return false;
         }
      }
      prime_count = primes->size();
      l.unlock();
      //
      // We're done on this segment, signal the next thread:
      //
      masks.reset();
      *current_max_processed_value = start_offset + max_points;
      start_offset += stride;
   }
   return true;
}

template<class ExecutionPolicy, class Integer, class Container, class OutputIterator>
void prime_sieve_imp(ExecutionPolicy&& policy, Integer upper_bound, Container& primes, OutputIterator out, bool output_to_container)
{
   if (upper_bound <= 2)
   {
      return;
   }

   bitmask_type sieve((upper_bound <= linear_sieve_limit<Integer> ? upper_bound : linear_sieve_limit<Integer>) / 2);

   if (output_to_container && (upper_bound > linear_sieve_limit<Integer>))
   {
      if (!jm::detail::linear_sieve_classical<Integer>(sieve, dual_output_iterator<OutputIterator, Container>(out, primes)))
         return;
   }
   else
   {
      if (!jm::detail::linear_sieve_classical<Integer>(sieve, out))
         return;
   }

   if (upper_bound <= linear_sieve_limit<Integer>)
      return;

   if constexpr (std::is_same_v<std::remove_cv_t<std::remove_reference_t<ExecutionPolicy>>, std::execution::sequenced_policy>
#if __cpp_lib_execution > 201900
      || std::is_same_v<std::remove_cv_t<std::remove_reference_t<ExecutionPolicy>>, std::execution::unsequenced_policy>
#endif
      )
   {
      for (Integer offset = linear_sieve_limit<Integer>; offset < upper_bound; offset += linear_sieve_limit<Integer>)
      {
         if (!linear_sieve_classical_segment(primes, sieve, offset, (std::min)(linear_sieve_limit<Integer>, upper_bound - offset), out, output_to_container))
            return;
      }
   }
   else
   {
      unsigned hardware_concurrency = std::thread::hardware_concurrency();
      if ((hardware_concurrency == 0) || (upper_bound <= linear_sieve_limit<Integer> * 2))
      {
         //
         // No point in using threads as there's only one more segment to process:
         //
         linear_sieve_classical_segment(primes, sieve, linear_sieve_limit<Integer>, upper_bound - linear_sieve_limit<Integer>, out, output_to_container);
      }
      else
      {
         unsigned n_threads = (upper_bound - linear_sieve_limit<Integer>) / linear_sieve_limit<Integer> +((upper_bound - linear_sieve_limit<Integer>) % linear_sieve_limit<Integer> ? 1 : 0);
         n_threads = (std::min)(n_threads, hardware_concurrency / 2);

         std::atomic<Integer> current_max_processed_value = linear_sieve_limit<Integer>;
         std::vector<std::unique_ptr<std::thread>> threads(n_threads);
         std::mutex mutex;

         for (unsigned i = 0; i < n_threads; ++i)
            threads[i].reset(new std::thread(linear_sieve_classical_segment_threaded<Container, Integer, OutputIterator>,
               &current_max_processed_value, &mutex, &primes, linear_sieve_limit<Integer> * (i + 1), upper_bound,
               linear_sieve_limit<Integer> * n_threads, out, output_to_container));

         for (unsigned i = 0; i < n_threads; ++i)
            threads[i]->join();
      }
   }
}

template <class OutputIterator>
struct counted_output_iterator
{
   std::shared_ptr<std::size_t> target;
   OutputIterator out;
   counted_output_iterator(std::size_t n, OutputIterator o) : target(new std::size_t()), out(o) 
   {
      *target = n;
   }

   struct counted_output_iterator_proxy
   {
      counted_output_iterator& out;
      counted_output_iterator_proxy(counted_output_iterator& o) : out(o) {}
      counted_output_iterator& operator*()const { return out; }
   };

   counted_output_iterator_proxy operator++(int){ return counted_output_iterator_proxy(*this); }
   counted_output_iterator& operator++() { return *this; }
   counted_output_iterator& operator*() { return *this; }

   template <class T>
   counted_output_iterator& operator=(const T& value) 
   { 
      *out++ = value; 
      --*target;
      return *this;
   }
};

template <class OutputIterator>
constexpr bool has_output_iterator_terminated(const counted_output_iterator<OutputIterator>& out)
{
   return *(out.target) == 0;
}

template <class Integer, class Container>
struct prime_factors_output_iterator
{
   Integer number;
   Container& out;

   prime_factors_output_iterator(Integer n, Container& o) : number(n), out(o) {}

   prime_factors_output_iterator& operator = (Integer prime)
   {
      std::size_t count = 0;
      while (number % prime == 0)
      {
         ++count;
         number /= prime;
         out.push_back(std::make_pair(prime, count));
      }
      return *this;
   }
   prime_factors_output_iterator& operator*() { return *this; }
   prime_factors_output_iterator& operator++() { return *this; }
   prime_factors_output_iterator& operator++(int) { return *this; }

   bool complete()const { return number == 1; }
};

template <class Integer, class Container>
constexpr bool has_output_iterator_terminated(const prime_factors_output_iterator<Integer, Container>& out)
{
   return out.complete();
}


} // End namespace detail

template<class Integer, class Container>
void prime_reserve(Integer upper_bound, Container &prime_container)
{
   typedef typename Container::size_type size_type;
   //
   // How big does the container have to be to hold all the primes in the range [2, upper_bound]?
   //
   // The prime number theorem (https://en.wikipedia.org/wiki/Prime_number_theorem) gives an asymptotic
   // estimate of upper_bound / log(upper_bound), but this only really holds as upper_bound -> infinity.
   //
   // Non asymptotic limits follow:
   //
   double limit;
   if (upper_bound >= 60184)
      // Dusart, Pierre (2010). "Estimates of Some Functions Over Primes without R.H"
      limit = upper_bound / (std::log(static_cast<double>(upper_bound)) - 1.1);
   else
   {
      // There are other loose relations available, but it's easy enough to search
      // all x in [3, 60184] and note how many primes there are and by how much
      // it differs from x / log(x).  The fudge factor of 1.25 is sufficient
      // always provide enough space:
      //
      limit = 1.25 * upper_bound / std::log(static_cast<double>(upper_bound));
   }
   prime_container.reserve(static_cast<size_type>(limit));
}

template<class ExecutionPolicy, class Integer, class Container>
inline std::enable_if_t<!boost::has_post_increment<Container>::value> prime_sieve(ExecutionPolicy&& policy, Integer upper_bound, Container &primes)
{
   typedef typename Container::value_type integer_type;
   prime_reserve(upper_bound, primes);
   return detail::prime_sieve_imp(policy, static_cast<integer_type>(upper_bound), primes, std::back_inserter(primes), false);
}

template<class ExecutionPolicy, class Integer, class OutputIterator>
inline std::enable_if_t<boost::has_post_increment<OutputIterator>::value> prime_sieve(ExecutionPolicy&& policy, Integer upper_bound, OutputIterator out)
{
   std::vector<Integer> primes;
   if (upper_bound > detail::linear_sieve_limit<Integer>)
      prime_reserve(upper_bound, primes);
   return detail::prime_sieve_imp(policy, upper_bound, primes, out, true);
}

template<class Integer, class Container>
inline std::enable_if_t<!boost::has_post_increment<Container>::value> prime_sieve(Integer upper_bound, Container &primes)
{
    prime_sieve(std::execution::seq, upper_bound, primes);
}
template<class Integer, class OutputIterator>
inline std::enable_if_t<boost::has_post_increment<OutputIterator>::value> prime_sieve(Integer upper_bound, OutputIterator primes)
{
    prime_sieve(std::execution::seq, upper_bound, primes);
}

template<class ExecutionPolicy, class Integer, class Container>
inline std::enable_if_t<!boost::has_post_increment<Container>::value> prime_sieve_n(ExecutionPolicy&& policy, Integer n, Container &primes)
{
   typedef typename Container::value_type integer_type;
   prime_reserve(n, primes);
   return detail::prime_sieve_imp(policy, (std::numeric_limits<Integer>::max)(), primes, detail::counted_output_iterator<std::back_insert_iterator<Container> >(n, std::back_inserter(primes)), false);
}

template<class ExecutionPolicy, class Integer, class OutputIterator>
inline std::enable_if_t<boost::has_post_increment<OutputIterator>::value> prime_sieve_n(ExecutionPolicy&& policy, Integer n, OutputIterator out)
{
   std::vector<Integer> primes;
   return detail::prime_sieve_imp(policy, (std::numeric_limits<Integer>::max)(), primes, detail::counted_output_iterator<OutputIterator>(n, out), true);
}

template<class Integer, class Container>
inline std::enable_if_t<!boost::has_post_increment<Container>::value> prime_sieve_n(Integer upper_bound, Container &primes)
{
    prime_sieve_n(std::execution::seq, upper_bound, primes);
}
template<class Integer, class OutputIterator>
inline std::enable_if_t<boost::has_post_increment<OutputIterator>::value> prime_sieve_n(Integer upper_bound, OutputIterator primes)
{
    prime_sieve_n(std::execution::seq, upper_bound, primes);
}

template<class ExecutionPolicy, class Integer>
inline std::vector<std::pair<Integer, std::size_t> > factorize(ExecutionPolicy&& policy, Integer n)
{
   std::vector<Integer> primes;
   std::vector<std::pair<Integer, std::size_t> > result;
   detail::prime_sieve_imp(policy, n / 2, primes, detail::prime_factors_output_iterator(n, result), true);
   return result;
}
template<class Integer>
inline std::vector<std::pair<Integer, std::size_t> > factorize(Integer n)
{
   return factorize(std::execution::seq, n);
}

}  // namespace boost::math


#endif //BOOST_MATH_SPECIAL_FUNCTIONS_PRIME_SIEVE_HPP
