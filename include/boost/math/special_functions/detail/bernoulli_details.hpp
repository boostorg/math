///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 John Maddock
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_BERNOULLI_DETAIL_HPP
#define BOOST_MATH_BERNOULLI_DETAIL_HPP

#include <boost/config.hpp>
#include <boost/detail/lightweight_mutex.hpp>

#ifdef BOOST_HAS_THREADS

#ifndef BOOST_NO_CXX11_HDR_ATOMIC
#  include <atomic>
#  define BOOST_MATH_ATOMIC_NS std
#if ATOMIC_INT_LOCK_FREE == 2
typedef std::atomic<int> atomic_counter_type;
#elif ATOMIC_SHORT_LOCK_FREE == 2
typedef std::atomic<short> atomic_counter_type;
#elif ATOMIC_LONG_LOCK_FREE == 2
typedef std::atomic<long> atomic_counter_type;
#elif ATOMIC_LLONG_LOCK_FREE == 2
typedef std::atomic<long long> atomic_counter_type;
#else
#  define BOOST_MATH_NO_ATOMIC_INT
#endif

#else // BOOST_NO_CXX11_HDR_ATOMIC
//
// We need Boost.Atomic, but on any platform that supports auto-linking we do
// not need to link against a separate library:
//
#define BOOST_ATOMIC_NO_LIB
#  include <boost/atomic.hpp>
#  define BOOST_MATH_ATOMIC_NS boost

namespace boost{ namespace math{ namespace detail{

//
// We need a type to use as an atomic counter:
//
#if BOOST_ATOMIC_INT_LOCK_FREE == 2
typedef boost::atomic<int> atomic_counter_type;
#elif BOOST_ATOMIC_SHORT_LOCK_FREE == 2
typedef boost::atomic<short> atomic_counter_type;
#elif BOOST_ATOMIC_LONG_LOCK_FREE == 2
typedef boost::atomic<long> atomic_counter_type;
#elif BOOST_ATOMIC_LLONG_LOCK_FREE == 2
typedef boost::atomic<long long> atomic_counter_type;
#else
#  define BOOST_MATH_NO_ATOMIC_INT
#endif

}}} // namespaces

#endif  // BOOST_NO_CXX11_HDR_ATOMIC

#endif // BOOST_HAS_THREADS

namespace boost{ namespace math{ namespace detail{

//
// We need to know the approximate value of /n/ which will
// cause bernoulli_b2n<T>(n) to return infinity - this allows
// us to elude a great deal of runtime checking for values below
// n, and only perform the full overflow checks when we know that we're
// getting close to the point where our calculations will overflow.
//
template <class T>
struct bernoulli_overflow_variant
{
   static const unsigned value = 
      (std::numeric_limits<T>::max_exponent == 128)
      && (std::numeric_limits<T>::radix == 2) ? 1 :
      (
         (std::numeric_limits<T>::max_exponent == 1024)
         && (std::numeric_limits<T>::radix == 2) ? 2 :
         (
            (std::numeric_limits<T>::max_exponent == 16384)
            && (std::numeric_limits<T>::radix == 2) ? 3 : 0
         )
      );
};

template<class T>
inline bool bernouli_impl_index_does_overflow(std::size_t, const mpl::int_<0>&)
{
   return false;  // We have no idea, perform a runtime check
}

// There are certain cases for which the index n, when trying
// to compute Bn, is known to overflow. In particular, when
// using 32-bit float and 64-bit double (IEEE 754 conformant),
// overflow will occur if the index exceeds the amount allowed
// in the tables of Bn.
template<class T>
inline bool bernouli_impl_index_does_overflow(std::size_t n, const mpl::int_<1>&)
{
   // This corresponds to 4-byte float, IEEE 745 conformant.
   return n >= max_bernoulli_index<1>::value * 2;
}
template<class T>
inline bool bernouli_impl_index_does_overflow(std::size_t n, const mpl::int_<2>&)
{
   // This corresponds to 8-byte float, IEEE 745 conformant.
   return n >= max_bernoulli_index<2>::value * 2;
}
template<class T>
inline bool bernouli_impl_index_does_overflow(std::size_t n, const mpl::int_<3>&)
{
   // This corresponds to 16-byte float, IEEE 745 conformant.
   return n >= max_bernoulli_index<3>::value * 2;
}

template<class T>
inline bool bernouli_impl_index_does_overflow(std::size_t n)
{
   typedef mpl::int_<bernoulli_overflow_variant<T>::value> tag_type;
   return bernouli_impl_index_does_overflow<T>(n, tag_type());
}

template<class T>
bool bernouli_impl_index_might_overflow(std::size_t n)
{
   if(bernouli_impl_index_does_overflow<T>(n))
   {
      // If the index *does* overflow, then it also *might* overflow.
      return true;
   }

   // Here, we use an asymptotic expansion of |Bn| from Luschny
   // to estimate if a given index n for Bn *might* overflow.
   static const T log_of_four_pi = log(boost::math::constants::two_pi<T>() * 2);
   static const T two_pi_e       = boost::math::constants::two_pi<T>() * boost::math::constants::e<T>();
   const float nf = static_cast<float>(n);
   const T nx (nf);
   const T nx2(nx * nx);
   const T n_log_term     = (nx + boost::math::constants::half<T>()) * log(nx / two_pi_e);

   const T approximate_log_of_bn =   log_of_four_pi
      + n_log_term
      + boost::math::constants::half<T>()
      + (T(1) / (nx * 12))
      - (T(1) / (nx2 * nx * 360))
      + (T(1) / (nx2 * nx2 * nx * 1260))
      + (n * boost::math::constants::ln_two<T>())
      - log(nx)
      + log(ldexp(T(1), (int)n) - 1);

   return approximate_log_of_bn * 1.1 > boost::math::tools::log_max_value<T>();
}

template<class T>
std::size_t possible_overflow_index()
{
   // We use binary search to determine a good approximation for an index that might overflow.
   // This code is called ONCE ONLY for each T at program startup.

   std::size_t upper_limit = 100000;
   std::size_t lower_limit = 1;

   if(bernouli_impl_index_might_overflow<T>(upper_limit * 2) == 0)
   {
      return upper_limit;
   }

   while(upper_limit > (lower_limit + 4))
   {
      const int mid = (upper_limit + lower_limit) / 2;

      if(bernouli_impl_index_might_overflow<T>(mid * 2) == 0)
      {
         lower_limit = mid;
      }
      else
      {
         upper_limit = mid;
      }
   }

   return lower_limit;
}
//
// The tangent numbers grow larger much more rapidly than the Bernoulli numbers do....
// so to compute the Bernoulli numbers from the tangent numbers, we need to avoid spurious
// overflow in the calculation, we can do this by scaling all the tangent number by some scale factor:
//
template <class T>
inline typename enable_if_c<std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::radix == 2), T>::type tangent_scale_factor()
{
   return ldexp(T(1), std::numeric_limits<T>::min_exponent + 5);
}
template <class T>
inline typename disable_if_c<std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::radix == 2), T>::type tangent_scale_factor()
{
   return tools::min_value<T>() * 16;
}
//
// Initializer: ensure all our constants are initialized prior to the first call of main:
//
template <class T, class Policy>
struct bernoulli_initializer
{
   struct init
   {
      init()
      {
         //
         // We call twice, once to initialize our static table, and once to
         // initialize our dymanic table:
         //
         boost::math::bernoulli_b2n<T>(2, Policy());
         boost::math::bernoulli_b2n<T>(max_bernoulli_b2n<T>::value + 1, Policy());
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
const typename bernoulli_initializer<T, Policy>::init bernoulli_initializer<T, Policy>::initializer;

//
// We need something to act as a cache for our calculated Bernoulli numbers.  In order to
// ensure both fast access and thread safety, we need a stable table which may be extended
// in size, but which never reallocates: that way values already calculated may be accessed
// concurrently with another thread extending the table with new values.
//
// Very very simple vector class that will never allocate more than once, we could use
// boost::container::static_vector here, but that allocates on the stack, which may well
// cause issues for the amount of memory we want in the extreme case...
//
template <class T>
struct fixed_vector : private std::allocator<T>
{
   typedef unsigned size_type;
   typedef T* iterator;
   typedef const T* const_iterator;
   fixed_vector() : m_used(0)
   { 
      std::size_t overflow_limit = 100 + possible_overflow_index<T>();
      m_capacity = (std::min)(overflow_limit, static_cast<std::size_t>(100000u));
      m_data = this->allocate(m_capacity); 
   }
   ~fixed_vector()
   {
      for(unsigned i = 0; i < m_used; ++i)
         this->destroy(&m_data[i]);
      this->deallocate(m_data, m_capacity);
   }
   T& operator[](unsigned n) { BOOST_ASSERT(n < m_used); return m_data[n]; }
   const T& operator[](unsigned n)const { BOOST_ASSERT(n < m_used); return m_data[n]; }
   unsigned size()const { return m_used; }
   unsigned size() { return m_used; }
   void resize(unsigned n, const T& val)
   {
      if(n > m_capacity)
         throw std::runtime_error("Exhausted storage for Bernoulli numbers.");
      for(unsigned i = m_used; i < n; ++i)
         new (m_data + i) T(val);
      m_used = n;
   }
   void resize(unsigned n) { resize(n, T()); }
   T* begin() { return m_data; }
   T* end() { return m_data + m_used; }
   T* begin()const { return m_data; }
   T* end()const { return m_data + m_used; }
private:
   T* m_data;
   unsigned m_used, m_capacity;
};

}}}

#endif // BOOST_MATH_BERNOULLI_DETAIL_HPP
