//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_SOLVE_HPP
#define BOOST_MATH_TOOLS_SOLVE_HPP

#ifdef _MSC_VER
#pragma once
#endif

#include <boost/config.hpp>
#include <boost/math/tools/assert.hpp>
#include <boost/math/tools/detail/linear_algebra/lu.hpp>

#include <cstddef>
#include <limits>
#include <vector>

namespace boost{ namespace math{ namespace tools{

// A minimal dense row-major matrix: all that solve() and remez.hpp need.
template <class T>
class matrix
{
public:
   matrix() : m_rows(0), m_cols(0) {}
   matrix(std::size_t rows, std::size_t cols) : m_rows(rows), m_cols(cols), m_data(rows * cols) {}
   T& operator()(std::size_t i, std::size_t j) { return m_data[i * m_cols + j]; }
   const T& operator()(std::size_t i, std::size_t j) const { return m_data[i * m_cols + j]; }
   std::size_t size1() const { return m_rows; }
   std::size_t size2() const { return m_cols; }
   const std::vector<T>& data() const { return m_data; }
private:
   std::size_t m_rows, m_cols;
   std::vector<T> m_data;
};

template <class T>
matrix<T> prod(const matrix<T>& A, const matrix<T>& B)
{
   BOOST_MATH_ASSERT(A.size2() == B.size1());
   matrix<T> result(A.size1(), B.size2());
   for(std::size_t i = 0; i < A.size1(); ++i)
      for(std::size_t j = 0; j < B.size2(); ++j)
      {
         T sum = 0;
         for(std::size_t k = 0; k < A.size2(); ++k)
            sum += A(i, k) * B(k, j);
         result(i, j) = sum;
      }
   return result;
}

template <class T>
std::vector<T> prod(const matrix<T>& A, const std::vector<T>& x)
{
   BOOST_MATH_ASSERT(A.size2() == x.size());
   std::vector<T> result(A.size1());
   for(std::size_t i = 0; i < A.size1(); ++i)
   {
      T sum = 0;
      for(std::size_t k = 0; k < A.size2(); ++k)
         sum += A(i, k) * x[k];
      result[i] = sum;
   }
   return result;
}

//
// Find x such that Ax = b
//
template <class T>
std::vector<T> solve(
          const matrix<T>& A_,
          const std::vector<T>& b_)
{
   BOOST_MATH_ASSERT(A_.size1() == b_.size() && A_.size2() == b_.size());

   const std::size_t n = b_.size();
   std::vector<T> A(A_.data()), b(b_), colmax(n);
   std::vector<std::size_t> piv(n);
   // The systems of the Remez method are close to singular, so only a zero pivot is refused:
   if(!detail::lu_factor(A, piv, colmax, n, T(0)))
      return std::vector<T>(n, std::numeric_limits<T>::quiet_NaN());
   detail::lu_solve(A, piv, n, b);
   //
   // iterate to reduce error:
   //
   std::vector<T> delta = prod(A_, b);
   for(std::size_t i = 0; i < delta.size(); ++i)
      delta[i] -= b_[i];
   detail::lu_solve(A, piv, n, delta);
   for(std::size_t i = 0; i < delta.size(); ++i)
      b[i] -= delta[i];

   return b;
}

}}} // namespaces

#endif // BOOST_MATH_TOOLS_SOLVE_HPP


