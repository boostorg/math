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
#include <boost/atomic.hpp>
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
T bernouli_impl_index_might_overflow_function(const T& nx)
{
   BOOST_MATH_STD_USING

   // Luschny LogB3(n) formula.

   const T nx2(nx * nx);

   const T t_one(1);

   const T approximate_log_of_bernoulli_bn
      =   ((boost::math::constants::half<T>() + nx) * log(nx))
        + ((boost::math::constants::half<T>() - nx) * log(boost::math::constants::pi<T>()))
        + (((T(3) / 2) - nx) * boost::math::constants::ln_two<T>())
        + ((nx * (T(2) - (nx2 * 7) * (t_one + ((nx2 * 30) * ((nx2 * 12) - t_one))))) / (((nx2 * nx2) * nx2) * 2520));

   return (approximate_log_of_bernoulli_bn - boost::math::tools::log_max_value<T>());
}

template<class T>
T bernouli_impl_index_might_overflow_tolerance(const T& min_value, const T& max_value)
{
   BOOST_MATH_STD_USING

   return abs(max_value - min_value) < T(0.9F);
}

template<class T>
inline std::size_t bernouli_impl_index_might_overflow()
{
   BOOST_MATH_STD_USING

   // Find the index of B2n that will overflow the numeric type.
   // Use Luschny's order-3 approximation of the logarithm of
   // Bernoulli numbers in combination with a root-finding tool
   // from boost::math::tools for this.

   boost::uintmax_t number_of_iterations = 128U;

   const std::pair<T, T> index_might_overflow_root =
      boost::math::tools::toms748_solve(bernouli_impl_index_might_overflow_function<T>,
                                        T(20),
                                        T((std::numeric_limits<std::size_t>::max)()),
                                        bernouli_impl_index_might_overflow_tolerance<T>,
                                        number_of_iterations);

   // Compute the index that might overflow by averaging the bracketed
   // values in the root pair.
   const T index_might_overflow_value = (index_might_overflow_root.first + index_might_overflow_root.second) / 2;

   // Determine a conservative estimate of the size_t integral value
   // of the index that might. Do this by evaluating the binary exponent
   // of the index that might overflow expressed as a floating-point number.
   int p2;
   frexp(index_might_overflow_value, &p2);

   // Armed with the binary exponent, we left shift by an amount
   // one less than this binary exponent such that the value of
   // the index that might overflow is guessed conservatively.

   const std::size_t index_might_overflow_conservative_value = static_cast<std::size_t>(static_cast<std::size_t>(1U) << (p2 - 1));

   // Example: A floating-point numeric type with a 32-bit signed
   // exponent overflows at an index of approximately 11,513,000.
   // The result of the conservative value for the overflow index
   // is approximately 8,389,000.

   return index_might_overflow_conservative_value;
}

template<class T>
std::size_t possible_overflow_index()
{
   // We use binary search to determine a good approximation for
   // an index that might overflow.
   // This code is called ONCE ONLY for each T at program startup.

   const std::size_t the_possible_overflow_index = bernouli_impl_index_might_overflow<T>();

   return the_possible_overflow_index;
}
//
// The tangent numbers grow larger much more rapidly than the Bernoulli numbers do....
// so to compute the Bernoulli numbers from the tangent numbers, we need to avoid spurious
// overflow in the calculation, we can do this by scaling all the tangent number by some scale factor:
//
template <class T>
inline typename enable_if_c<std::numeric_limits<T>::is_specialized && (std::numeric_limits<T>::radix == 2), T>::type tangent_scale_factor()
{
   BOOST_MATH_STD_USING
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
         try{
            boost::math::bernoulli_b2n<T>(max_bernoulli_b2n<T>::value + 1, Policy());
         } catch(const std::overflow_error&){}
         boost::math::tangent_t2n<T>(2, Policy());
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
      std::size_t overflow_limit = 100 + 3 * possible_overflow_index<T>();
      m_capacity = static_cast<size_type>((std::min)(overflow_limit, static_cast<std::size_t>(100000u)));
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

template <class T, class Policy>
class bernoulli_numbers_cache
{
public:
   bernoulli_numbers_cache() : m_overflow_limit((std::numeric_limits<std::size_t>::max)())
#if defined(BOOST_HAS_THREADS) && !defined(BOOST_MATH_NO_ATOMIC_INT)
      , m_counter(0)
#endif
   {}

   typedef fixed_vector<T> container_type;

   void tangent(std::size_t m)
   {
      static const std::size_t min_overflow_index = possible_overflow_index<T>();
      tn.resize(static_cast<container_type::size_type>(m), T(0U));

      BOOST_MATH_INSTRUMENT_VARIABLE(min_overflow_index);

      std::size_t prev_size = m_intermediates.size();
      m_intermediates.resize(m, T(0U));

      if(prev_size == 0)
      {
         m_intermediates[1] = tangent_scale_factor<T>() /*T(1U)*/;
         tn[0U] = T(0U);
         tn[1U] = tangent_scale_factor<T>()/* T(1U)*/;
         BOOST_MATH_INSTRUMENT_VARIABLE(tn[0]);
         BOOST_MATH_INSTRUMENT_VARIABLE(tn[1]);
      }

      for(std::size_t i = std::max<size_t>(2, prev_size); i < m; i++)
      {
         bool overflow_check = false;
         if(i >= min_overflow_index && (boost::math::tools::max_value<T>() / (i-1) < m_intermediates[1]) )
         {
            std::fill(tn.begin() + i, tn.end(), boost::math::tools::max_value<T>());
            break;
         }
         m_intermediates[1] = m_intermediates[1] * (i-1);
         for(std::size_t j = 2; j <= i; j++)
         {
            overflow_check =
                  (i >= min_overflow_index) && (
                  (boost::math::tools::max_value<T>() / (i - j) < m_intermediates[j])
                  || (boost::math::tools::max_value<T>() / (i - j + 2) < m_intermediates[j-1])
                  || (boost::math::tools::max_value<T>() - m_intermediates[j] * (i - j) < m_intermediates[j-1] * (i - j + 2))
                  || ((boost::math::isinf)(m_intermediates[j]))
                );

            if(overflow_check)
            {
               std::fill(tn.begin() + i, tn.end(), boost::math::tools::max_value<T>());
               break;
            }
            m_intermediates[j] = m_intermediates[j] * (i - j) + m_intermediates[j-1] * (i - j + 2);
         }
         if(overflow_check)
            break; // already filled the tn...
         tn[static_cast<container_type::size_type>(i)] = m_intermediates[i];
         BOOST_MATH_INSTRUMENT_VARIABLE(i);
         BOOST_MATH_INSTRUMENT_VARIABLE(tn[i]);
      }
   }

   void tangent_numbers_series(const std::size_t m)
   {
      BOOST_MATH_STD_USING
      static const std::size_t min_overflow_index = possible_overflow_index<T>();

      typename container_type::size_type old_size = bn.size();

      tangent(m);
      bn.resize(static_cast<container_type::size_type>(m));

      if(!old_size)
      {
         bn[0] = 1;
         old_size = 1;
      }

      T power_two(ldexp(T(1), static_cast<int>(2 * old_size)));

      for(std::size_t i = old_size; i < m; i++)
      {
         T b(static_cast<T>(i * 2));
         //
         // Not only do we need to take care to avoid spurious over/under flow in
         // the calculation, but we also need to avoid overflow altogether in case
         // we're calculating with a type where "bad things" happen in that case:
         //
         b  = b / (power_two * tangent_scale_factor<T>());
         b /= (power_two - 1);
         bool overflow_check = (i >= min_overflow_index) && (tools::max_value<T>() / tn[static_cast<container_type::size_type>(i)] < b);
         if(overflow_check)
         {
            m_overflow_limit = i;
            while(i < m)
            {
               b = std::numeric_limits<T>::has_infinity ? std::numeric_limits<T>::infinity() : tools::max_value<T>();
               bn[static_cast<container_type::size_type>(i)] = ((i % 2U) ? b : -b);
               ++i;
            }
            break;
         }
         else
         {
            b *= tn[static_cast<container_type::size_type>(i)];
         }

         power_two = ldexp(power_two, 2);

         const bool b_neg = i % 2 == 0;

         bn[static_cast<container_type::size_type>(i)] = ((!b_neg) ? b : -b);
      }
   }

   template <class OutputIterator>
   OutputIterator copy_bernoulli_numbers(OutputIterator out, std::size_t start, std::size_t n, const Policy& pol)
   {
      //
      // There are basically 3 thread safety options:
      //
      // 1) There are no threads (BOOST_HAS_THREADS is not defined).
      // 2) There are threads, but we do not have a true atomic integer type, 
      //    in this case we just use a mutex to guard against race conditions.
      // 3) There are threads, and we have an atomic integer: in this case we can
      //    use the double-checked locking pattern to avoid thread synchronisation
      //    when accessing values already in the cache.
      //
   #if !defined(BOOST_HAS_THREADS)
      //
      // Single threaded code, very simple:
      //
      if(start + n >= bn.size())
      {
         std::size_t new_size = (std::max)((std::max)(start + n, std::size_t(bn.size() + 20)), std::size_t(50));
         tangent_numbers_series(new_size);
      }

      for(std::size_t i = (std::max)(max_bernoulli_b2n<T>::value + 1, start); i < start + n; ++i)
      {
         *out = (i >= m_overflow_limit) ? policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol) : bn[i];
         ++out;
      }
   #elif defined(BOOST_MATH_NO_ATOMIC_INT)
      //
      // We need to grab a mutex every time we get here, for both readers and writers:
      //
      boost::detail::lightweight_mutex::scoped_lock l(m_mutex);
      if(start + n >= bn.size())
      {
         std::size_t new_size = (std::max)((std::max)(start + n, std::size_t(bn.size() + 20)), std::size_t(50));
         tangent_numbers_series(new_size);
      }

      for(std::size_t i = (std::max)(max_bernoulli_b2n<T>::value + 1, start); i < start + n; ++i)
      {
         *out = (i >= m_overflow_limit) ? policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol) : bn[i];
         ++out;
      }

   #else
      //
      // Double-checked locking pattern, lets us access cached already cached values
      // without locking:
      //
      // Get the counter and see if we need to calculate more constants:
      //
      if(static_cast<std::size_t>(m_counter.load(BOOST_MATH_ATOMIC_NS::memory_order_consume)) < start + n)
      {
         boost::detail::lightweight_mutex::scoped_lock l(m_mutex);

         if(static_cast<std::size_t>(m_counter.load(BOOST_MATH_ATOMIC_NS::memory_order_consume)) < start + n)
         {
            if(start + n >= bn.size())
            {
               std::size_t new_size = (std::max)((std::max)(start + n, std::size_t(bn.size() + 20)), std::size_t(50));
               tangent_numbers_series(new_size);
            }
            m_counter.store(bn.size(), BOOST_MATH_ATOMIC_NS::memory_order_release);
         }
      }

      for(std::size_t i = (std::max)(static_cast<std::size_t>(max_bernoulli_b2n<T>::value + 1), start); i < start + n; ++i)
      {
         *out = (i >= m_overflow_limit) ? policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol) : bn[static_cast<container_type::size_type>(i)];
         ++out;
      }

   #endif
      return out;
   }

   template <class OutputIterator>
   OutputIterator copy_tangent_numbers(OutputIterator out, std::size_t start, std::size_t n, const Policy& pol)
   {
      //
      // There are basically 3 thread safety options:
      //
      // 1) There are no threads (BOOST_HAS_THREADS is not defined).
      // 2) There are threads, but we do not have a true atomic integer type, 
      //    in this case we just use a mutex to guard against race conditions.
      // 3) There are threads, and we have an atomic integer: in this case we can
      //    use the double-checked locking pattern to avoid thread synchronisation
      //    when accessing values already in the cache.
      //
   #if !defined(BOOST_HAS_THREADS)
      //
      // Single threaded code, very simple:
      //
      if(start + n >= bn.size())
      {
         std::size_t new_size = (std::max)((std::max)(start + n, std::size_t(bn.size() + 20)), std::size_t(50));
         tangent_numbers_series(new_size);
      }

      for(std::size_t i = start; i < start + n; ++i)
      {
         if(i >= m_overflow_limit)
            *out = policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol);
         else
         {
            if(tools::max_value<T>() * tangent_scale_factor<T>() < tn[i])
               *out = policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol);
            else
               *out = tn[i] / tangent_scale_factor<T>();
         }
         ++out;
      }
   #elif defined(BOOST_MATH_NO_ATOMIC_INT)
      //
      // We need to grab a mutex every time we get here, for both readers and writers:
      //
      boost::detail::lightweight_mutex::scoped_lock l(m_mutex);
      if(start + n >= bn.size())
      {
         std::size_t new_size = (std::max)((std::max)(start + n, std::size_t(bn.size() + 20)), std::size_t(50));
         tangent_numbers_series(new_size);
      }

      for(std::size_t i = start; i < start + n; ++i)
      {
         if(i >= m_overflow_limit)
            *out = policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol);
         else
         {
            if(tools::max_value<T>() * tangent_scale_factor<T>() < tn[i])
               *out = policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol);
            else
               *out = tn[i] / tangent_scale_factor<T>();
         }
         ++out;
      }

   #else
      //
      // Double-checked locking pattern, lets us access cached already cached values
      // without locking:
      //
      // Get the counter and see if we need to calculate more constants:
      //
      if(static_cast<std::size_t>(m_counter.load(BOOST_MATH_ATOMIC_NS::memory_order_consume)) < start + n)
      {
         boost::detail::lightweight_mutex::scoped_lock l(m_mutex);

         if(static_cast<std::size_t>(m_counter.load(BOOST_MATH_ATOMIC_NS::memory_order_consume)) < start + n)
         {
            if(start + n >= bn.size())
            {
               std::size_t new_size = (std::max)((std::max)(start + n, std::size_t(bn.size() + 20)), std::size_t(50));
               tangent_numbers_series(new_size);
            }
            m_counter.store(bn.size(), BOOST_MATH_ATOMIC_NS::memory_order_release);
         }
      }

      for(std::size_t i = start; i < start + n; ++i)
      {
         if(i >= m_overflow_limit)
            *out = policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol);
         else
         {
            if(tools::max_value<T>() * tangent_scale_factor<T>() < tn[static_cast<container_type::size_type>(i)])
               *out = policies::raise_overflow_error<T>("boost::math::bernoulli<%1%>(std::size_t)", 0, pol);
            else
               *out = tn[static_cast<container_type::size_type>(i)] / tangent_scale_factor<T>();
         }
         ++out;
      }

   #endif
      return out;
   }

private:
   //
   // The caches for Bernoulli and tangent numbers, once allocated,
   // these must NEVER EVER reallocate as it breaks our thread
   // safety guarentees:
   //
   fixed_vector<T> bn, tn;
   std::vector<T> m_intermediates;
   // The value at which we know overflow has already occured for the Bn:
   std::size_t m_overflow_limit;
#if !defined(BOOST_HAS_THREADS)
#elif defined(BOOST_MATH_NO_ATOMIC_INT)
   boost::detail::lightweight_mutex m_mutex;
#else
   boost::detail::lightweight_mutex m_mutex;
   atomic_counter_type m_counter;
#endif
};

template <class T, class Policy>
inline bernoulli_numbers_cache<T, Policy>& get_bernoulli_numbers_cache()
{
   //
   // Force this function to be called at program startup so all the static variables
   // get initailzed then (thread safety).
   //
   bernoulli_initializer<T, Policy>::force_instantiate();
   static bernoulli_numbers_cache<T, Policy> data;
   return data;
}

}}}

#endif // BOOST_MATH_BERNOULLI_DETAIL_HPP
