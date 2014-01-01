///////////////////////////////////////////////////////////////////////////////
//  Copyright 2013 John Maddock
//  Copyright 2013 Christopher Kormanyos
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FIXED_CAPACITY_HEAP_VECTOR_HPP
#define BOOST_MATH_FIXED_CAPACITY_HEAP_VECTOR_HPP

namespace boost { namespace math { namespace tools {

#include <iterator>
#include <memory>
#include <boost/assert.hpp>

//
// The fixed_capacity_heap_vector container acts as a cache for dynamically
// calculated numbers. In order to ensure both fast access and thread safety,
// this container may be extended in size, but never reallocates beyond
// its maximum capacity. The maximum capacity can be given as a template
// parameter, whereby a smaller value can be used in the class constructor.
// That way values already calculated may be accessed concurrently with
// another thread extending the table with new values.
//
// This is a very, very simple vector class that will never allocate
// more than once, we could use boost::container::static_vector here,
// but that allocates on the stack, which may well cause issues for
// the amount of memory we want in the extreme case.
//

template <class T, const unsigned maximum_capacity = 100000U>
struct fixed_capacity_heap_vector : private std::allocator<T>
{
   typedef T                                     value_type;
   typedef unsigned                              size_type;
   typedef value_type&                           reference;
   typedef const value_type&                     const_reference;
   typedef value_type*                           iterator;
   typedef const value_type*                     const_iterator;
   typedef std::reverse_iterator<iterator>       reverse_iterator;
   typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

   fixed_capacity_heap_vector(const size_type this_capacity)
      : m_capacity((std::min)((std::max)(size_type(100U), this_capacity), size_type(maximum_capacity))),
        m_used    (0U),
        m_data    (this->allocate(m_capacity)) { }

   ~fixed_capacity_heap_vector()
   {
      for(size_type i = 0U; i < m_used; ++i)
         this->destroy(&m_data[i]);

      this->deallocate(m_data, m_capacity);
   }

         reference operator[](size_type n)       { BOOST_ASSERT(n < m_used); return m_data[n]; }
   const_reference operator[](size_type n) const { BOOST_ASSERT(n < m_used); return m_data[n]; }

   reference at(const size_type i)             { BOOST_ASSERT(n < m_used); return m_data[n]; }
   const_reference at(const size_type i) const { BOOST_ASSERT(n < m_used); return m_data[n]; }

   reference front()             { return m_data[0U]; }
   const_reference front() const { return m_data[0U]; }

   reference back()             { return m_data[m_used - 1U]; }
   const_reference back() const { return m_data[m_used - 1U]; }

   size_type size() const { return m_used; }

   size_type max_size() const { return m_capacity; }

   bool empty() const { return (m_used == size_type(0U)); }

   void resize(size_type n, const value_type& val)
   {
      if(n > m_capacity)
         throw std::runtime_error("fixed_capacity_heap_vector is out of storage.");

      for(size_type i = m_used; i < n; ++i)
         new (m_data + i) T(val);

      m_used = n;
   }

   void resize(size_type n) { resize(n, T()); }

   iterator begin() { return m_data; }
   iterator end  () { return m_data + m_used; }

   const_iterator begin() const { return m_data; }
   const_iterator end  () const { return m_data + m_used; }

   const_iterator cbegin() const { return m_data; }
   const_iterator cend  () const { return m_data + m_used; }

   reverse_iterator rbegin() { return reverse_iterator(m_data + m_used); }
   reverse_iterator rend  () { return reverse_iterator(m_data); }

   const_reverse_iterator rbegin() const { return const_reverse_iterator(m_data + m_used); }
   const_reverse_iterator rend  () const { return const_reverse_iterator(m_data); }

   const_reverse_iterator crbegin() const { return const_reverse_iterator(m_data + m_used); }
   const_reverse_iterator crend  () const { return const_reverse_iterator(m_data); }

private:
   const size_type m_capacity;
   size_type m_used;
   iterator m_data;
};

}}} // namespace boost::math::tools

#endif // BOOST_MATH_FIXED_CAPACITY_HEAP_VECTOR_HPP
