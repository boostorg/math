//  (C) Copyright John Maddock 2006.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_POLYNOMIAL_HPP
#define BOOST_MATH_TOOLS_POLYNOMIAL_HPP

#include <boost/assert.hpp>
#include <boost/math/tools/rational.hpp>
#include <boost/math/tools/real_cast.hpp>

#include <vector>

namespace boost{ namespace math{ namespace tools{

template <class T>
class polynomial
{
public:
   // typedefs:
   typedef typename std::vector<T>::value_type value_type;
   typedef typename std::vector<T>::size_type size_type;

   // construct:
   polynomial(){}
   template <class U>
   polynomial(const U* data, unsigned order)
      : m_data(data, data + order + 1)
   {
   }
   template <class U>
   polynomial(const U& point)
   {
      m_data.push_back(point);
   }

   // copy:
   polynomial(const polynomial& p)
      : m_data(p.m_data) { }

   template <class U>
   polynomial(const polynomial<U>& p)
   {
      for(unsigned i = 0; i < p.size(); ++i)
      {
         m_data.push_back(boost::math::tools::real_cast<T>(p[i]));
      }
   }

   // access:
   size_type size()const { return m_data.size(); }
   size_type degree()const { return m_data.size() - 1; }
   value_type& operator[](size_type i)
   {
      return m_data[i];
   }
   const value_type& operator[](size_type i)const
   {
      return m_data[i];
   }
   T evaluate(T z)const
   {
      return boost::math::tools::evaluate_polynomial(&m_data[0], z, m_data.size());;
   }

   // operators:
   template <class U>
   polynomial& operator +=(const U& value)
   {
      if(m_data.size() == 0)
         m_data.push_back(value);
      else
      {
         m_data[0] += value;
      }
      return *this;
   }
   template <class U>
   polynomial& operator -=(const U& value)
   {
      if(m_data.size() == 0)
         m_data.push_back(-value);
      else
      {
         m_data[0] -= value;
      }
      return *this;
   }
   template <class U>
   polynomial& operator *=(const U& value)
   {
      for(size_type i = 0; i < m_data.size(); ++i)
         m_data[i] *= value;
      return *this;
   }
   template <class U>
   polynomial& operator +=(const polynomial<U>& value)
   {
      size_type s1 = (std::min)(m_data.size(), value.size());
      for(size_type i = 0; i < s1; ++i)
         m_data[i] += value[i];
      for(size_type i = s1; i < value.size(); ++i)
         m_data.push_back(value[i]);
      return *this;
   }
   template <class U>
   polynomial& operator -=(const polynomial<U>& value)
   {
      size_type s1 = (std::min)(m_data.size(), value.size());
      for(size_type i = 0; i < s1; ++i)
         m_data[i] -= value[i];
      for(size_type i = s1; i < value.size(); ++i)
         m_data.push_back(-value[i]);
      return *this;
   }
   template <class U>
   polynomial& operator *=(const polynomial<U>& value)
   {
      // TODO: FIXME: use O(N log(N)) algorithm!!!
      BOOST_ASSERT(value.size());
      polynomial base(*this);
      *this *= value[0];
      for(size_type i = 1; i < value.size(); ++i)
      {
         polynomial t(base);
         t *= value[i];
         size_type s = size() - i;
         for(size_type j = 0; j < s; ++j)
         {
            m_data[i+j] += t[j];
         }
         for(size_type j = s; j < t.size(); ++j)
            m_data.push_back(t[j]);
      }
      return *this;
   }

private:
   std::vector<T> m_data;
};

template <class T>
inline polynomial<T> operator + (const polynomial<T>& a, const polynomial<T>& b)
{
   polynomial<T> result(a);
   result += b;
   return result;
}

template <class T>
inline polynomial<T> operator - (const polynomial<T>& a, const polynomial<T>& b)
{
   polynomial<T> result(a);
   result -= b;
   return result;
}

template <class T>
inline polynomial<T> operator * (const polynomial<T>& a, const polynomial<T>& b)
{
   polynomial<T> result(a);
   result *= b;
   return result;
}

template <class T, class U>
inline polynomial<T> operator + (const polynomial<T>& a, const U& b)
{
   polynomial<T> result(a);
   result += b;
   return result;
}

template <class T, class U>
inline polynomial<T> operator - (const polynomial<T>& a, const U& b)
{
   polynomial<T> result(a);
   result -= b;
   return result;
}

template <class T, class U>
inline polynomial<T> operator * (const polynomial<T>& a, const U& b)
{
   polynomial<T> result(a);
   result *= b;
   return result;
}

template <class U, class T>
inline polynomial<T> operator + (const U& a, const polynomial<T>& b)
{
   polynomial<T> result(b);
   result += a;
   return result;
}

template <class U, class T>
inline polynomial<T> operator - (const U& a, const polynomial<T>& b)
{
   polynomial<T> result(a);
   result -= b;
   return result;
}

template <class U, class T>
inline polynomial<T> operator * (const U& a, const polynomial<T>& b)
{
   polynomial<T> result(b);
   result *= a;
   return result;
}

template <class charT, class traits, class T>
inline std::basic_ostream<charT, traits>& operator << (std::basic_ostream<charT, traits>& os, const polynomial<T>& poly)
{
   os << "{ ";
   for(unsigned i = 0; i < poly.size(); ++i)
   {
      if(i) os << ", ";
      os << poly[i];
   }
   os << " }";
   return os;
}

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_POLYNOMIAL_HPP


