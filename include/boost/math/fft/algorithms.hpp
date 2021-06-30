///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FFT_ALGORITHMS_HPP
  #define BOOST_MATH_FFT_ALGORITHMS_HPP

  #include <algorithm>
  #include <numeric>
  #include <cmath>
  #include <vector>
  #include <boost/math/constants/constants.hpp>
  #include <boost/math/fft/abstract_ring.hpp>
  


  namespace boost { namespace math {  namespace fft {

  namespace detail {
  
  template <typename integer>
  integer modulo(const integer a, const integer b)
  {
    return ( a%b + b ) % b;
  }
  
  template <typename integer>
  struct euclid_return
  {
    integer gcd, x, y;
  };

  template <typename integer>
  euclid_return<integer> extended_euclid(const integer a, const integer b)
  {
   static_assert(std::numeric_limits<integer>::is_signed,
                 "Signed integers are needed");
   integer x = 1, y = 0, g = a;
   for (integer x1 = 0, y1 = 1, g1 = b; g1;)
   {
     integer q = g / g1;
     x = x - x1 * q;
     std::swap(x, x1);
     y = y - y1 * q;
     std::swap(y, y1);
     g = g - g1 * q;
     std::swap(g, g1);
   }
   if(g<0)
   {
    g = -g;
    x = -x;
    y = -y;
   }
   return {g, x, y};
  }
  
  template <typename integer>
  integer gcd(const integer a, const integer b)
  {
    euclid_return<integer> ret = extended_euclid(a,b);
    return ret.gcd;
  }
  
  
  // TODO: use another power function
  template <class T>
  T power(const T& x, int n)
  {
    /*
        for our use case n should always be >0.
        However, if n==0 we would expect something
        like 1.
    */
    // TODO: decide which exceptions to throw
    if (n < 1)
        throw std::domain_error("power(x,n) expects n>=1.");

    bool identity = true;
    T r{}, aux{x};

    for (; n; n >>= 1)
    {
      if (n & 1)
      {
        r = identity ? aux : r * aux;
        identity = false;
      }
      aux *= aux;
    }
    return r;
  }
    
  template<class ComplexType>
  ComplexType complex_root_of_unity(long n,long p=1)
  /*
    Computes exp(-i 2 pi p/n)
  */
  {
    using real_value_type = typename ComplexType::value_type;
    p = modulo(p,n);
    
    if(p==0)
      return ComplexType(1,0);
    
    long g = gcd(p,n); 
    n/=g;
    p/=g;
    switch(n)
    {
      case 1:
        return ComplexType(1,0);
      case 2:
        return p==0 ? ComplexType(1,0) : ComplexType(-1,0);
      case 4:
        return p==0 ? ComplexType(1,0) : 
               p==1 ? ComplexType(0,-1) :
               p==2 ? ComplexType(-1,0) :
                      ComplexType(0,1) ;
    }
    using std::sin;
    using std::cos;
    real_value_type phase = -2*p*boost::math::constants::pi<real_value_type>()/n;
    return ComplexType(cos(phase),sin(phase));
  }
  template<class ComplexType>
  ComplexType complex_inverse_root_of_unity(long n,long p=1)
  /*
    Computes exp(i 2 pi p/n)
  */
  {
    return complex_root_of_unity<ComplexType>(n,-p);
  }
  inline bool is_power2(long x) { return x == (x & -x);}
  
  inline long upper_bound_power2(long x)
  {
    /*
      Returns the smallert power of two that is greater or equal to x.
    */
    long up=1;
    for(;up<x;up<<=1);
    return up;
  }
  
  inline long lower_bound_power2(long x)
  {
    /*
      Returns the biggest power of two that is smaller than or equal to x.
    */
    for(long lsb = x & -x; lsb!=x;)
    {
      x-=lsb;
      lsb = x&-x;
    }
    return x;
  }
    
  std::vector<long> prime_factorization(long n)
  {
    /*
        Naive O(sqrt(n)) prime factorization.
    */
    std::vector<long> F;
    for (long x = 2; x * x <= n;)
      if (n % x == 0)
      {
        F.push_back(x);
        n /= x;
      }
      else
      {
        ++x;
      }
    if (n > 1)
      F.push_back(n);
    return F;
  }
  
  
  template<class T>
  void dft_prime_bruteForce(const T* in_first, const T* in_last, T* out, const T w)
  /*
    assumptions: 
    - allocated memory in out is enough to hold distance(in_first,in_last) element,
  */
  {
    const long N = static_cast<long>(std::distance(in_first,in_last));
    if(N<=0)
      return;
    
    std::unique_ptr<T> mem;
    T* work_space = out;
    
    if(in_first == out)
    {
      mem.reset(new T[N]);
      work_space = mem.get();
    }
    
    work_space[0] = std::accumulate(in_first+1,in_last,in_first[0]);
    
    T wi=w;
    for(long i=1;i<N;++i, wi*=w)
    {
      T wij{wi};
      T sum{in_first[0]};
      for(long j=1;j<N; ++j, wij*=wi)
      {
        sum += in_first[j]*wij;
      }
      work_space[i] = sum;
    }
    
    if(out != work_space)
      std::copy(work_space,work_space+N,out);
  }
  
  template<class complex_value_type>
  void complex_dft_prime_bruteForce(
    const complex_value_type* in_first, 
    const complex_value_type* in_last, 
    complex_value_type* out, int sign)
  /*
    assumptions: 
    - allocated memory in out is enough to hold distance(in_first,in_last) element,
  */
  {
    const long N = static_cast<long>(std::distance(in_first,in_last));
    if(N<=0)
      return;
    
    std::unique_ptr<complex_value_type> mem;
    complex_value_type* work_space = out;
    
    if(in_first == out)
    {
      mem.reset(new complex_value_type[N]);
      work_space = mem.get();
    }
    
    work_space[0] = std::accumulate(in_first+1,in_last,in_first[0]);
    
    
    for(long i=1;i<N;++i)
    {
      complex_value_type sum{in_first[0]};
      for(long j=1;j<N; ++j)
      {
        sum += in_first[j] * complex_root_of_unity<complex_value_type>(N,i*j*sign);
      }
      work_space[i] = sum;
    }
    
    if(out != work_space)
      std::copy(work_space,work_space+N,out);
  }
  
  
  template <class T>
  void dft_composite(const T *in_first, const T *in_last, T* out, const T e)
  {
    /*
      Cooley-Tukey mapping, intrinsically out-of-place, Decimation in Time
      composite sizes.
    */
    const long n = static_cast<unsigned int>(std::distance(in_first,in_last));
    if(n <=0 )
      return;
    
    if (n == 1)
    {
        out[0]=in_first[0];
        return;
    }
    auto prime_factors = prime_factorization(n);
    
    // reorder input
    for (long i = 0; i < n; ++i)
    {
        long j = 0, k = i;
        for (auto p : prime_factors)
        {
            j = j * p + k % p;
            k /= p;
        }
        out[j] = in_first[i];
    }
    
    std::reverse(prime_factors.begin(), prime_factors.end());
    
    // butterfly pattern
    long len = 1;
    for (auto p : prime_factors)
    {
      long len_old = len;
      len *= p;
      T w_len = power(e, n / len);
      T w_p = power(e,n/p);
      
      std::vector<T> tmp(p);
      for (long i = 0; i < n; i += len)
      {
        for(long k=0;k<len_old;++k)
        {
          for(long j=0;j<p;++j)
            if(j==0 or k==0)
              tmp[j] = out[i + j*len_old +k ];
            else
              tmp[j] = out[i + j*len_old +k ] * power(w_len,k*j);
          
          dft_prime_bruteForce(tmp.data(),tmp.data()+p,tmp.data(),w_p);
          
          for(long j=0;j<p;++j)
            out[i+ j*len_old + k] = tmp[j];
        }
      }
    }
  }
  
  template <class ComplexType>
  void complex_dft_composite(const ComplexType *in_first, const ComplexType *in_last, ComplexType* out, int sign)
  {
    /*
      Cooley-Tukey mapping, intrinsically out-of-place, Decimation in Time
      composite sizes.
    */
    const long n = static_cast<long>(std::distance(in_first,in_last));
    if(n <=0 )
      return;
    
    if (n == 1)
    {
        out[0]=in_first[0];
        return;
    }
    auto prime_factors = prime_factorization(n);
    
    // reorder input
    for (long i = 0; i < n; ++i)
    {
        long j = 0, k = i;
        for (auto p : prime_factors)
        {
            j = j * p + k % p;
            k /= p;
        }
        out[j] = in_first[i];
    }
    
    std::reverse(prime_factors.begin(), prime_factors.end());
    
    // butterfly pattern
    long len = 1;
    for (auto p : prime_factors)
    {
      long len_old = len;
      len *= p;
      
      std::vector<ComplexType> tmp(p);
      for (long i = 0; i < n; i += len)
      {
        for(long k=0;k<len_old;++k)
        {
          for(long j=0;j<p;++j)
            if(j==0 || k==0)
              tmp[j] = out[i + j*len_old +k ];
            else
              tmp[j] = out[i + j*len_old +k ] * complex_root_of_unity<ComplexType>(len,k*j*sign);
          
          complex_dft_prime_bruteForce(tmp.data(),tmp.data()+p,tmp.data(),sign);
          
          for(long j=0;j<p;++j)
            out[i+ j*len_old + k] = tmp[j];
        }
      }
    }
  }
  
  template <class T>
  void dft_power2(const T *in_first, const T *in_last, T* out, const T e)
  {
    /*
      Cooley-Tukey mapping, in-place Decimation in Time 
    */
    const long ptrdiff = std::distance(in_first,in_last);
    if(ptrdiff <=0 )
      return;
    const long n = lower_bound_power2(ptrdiff);
    
    if(in_first!=out)
      std::copy(in_first,in_last,out);
    
    if (n == 1)
        return;

    // auto _1 = T{1};
    
    int nbits = 0;
    std::vector<T> e2{e};
    for (int m = n / 2; m > 0; m >>= 1, ++nbits)
      e2.push_back(e2.back() * e2.back());

    std::reverse(e2.begin(), e2.end());

    // Gold-Rader bit-reversal algorithm.
    for(int i=0,j=0;i<n-1;++i)
    { 
      if(i<j)
        std::swap(out[i],out[j]);
      for(int k=n>>1;!( (j^=k)&k );k>>=1);
    }
    
    
    for (int len = 2, k = 1; len <= n; len <<= 1, ++k)
    {
      for (int i = 0; i < n; i += len)
      {
        {
          int j=0;
          T* u = out + i + j, *v = out + i + j + len / 2;
          T Bu = *u, Bv = *v;
          *u = Bu + Bv;
          *v = Bu - Bv;
        }
        
        T ej = e2[k];
        for (int j = 1; j < len / 2; ++j)
        {
          T* u = out + i + j, *v = out + i + j + len / 2;
          T Bu = *u, Bv = *v * ej;
          *u = Bu + Bv;
          *v = Bu - Bv;
          ej *= e2[k];
        }
      }
    }
  }
  
  template<class complex_value_type>
  void complex_dft_power2(
    const complex_value_type *in_first, 
    const complex_value_type *in_last, 
    complex_value_type* out, int sign)
  {
    // Naive in-place complex DFT.
    const long ptrdiff = static_cast<long>(std::distance(in_first,in_last));
    if(ptrdiff <=0 )
      return;
    const long my_n = lower_bound_power2(ptrdiff);
    
    if(in_first!=out)
      std::copy(in_first,in_last,out);
    
    if (my_n == 1)
        return;
    
    if(in_first!=out)
      std::copy(in_first, in_last, out);

    // Recursive decimation in frequency.
    for(long m = my_n; m > 1; m /= 2)
    {
      long mh = m / 2;
      
      for(long j = 0; j < mh; ++j)
      {
        complex_value_type cs{complex_root_of_unity<complex_value_type>(m,j*sign)};

        for (long t1=j; t1 < j + my_n; t1 += m)
        {
          complex_value_type u = out[t1];
          complex_value_type v = out[t1 + mh];

          out[t1]      =  u + v;
          out[t1 + mh] = (u - v) * cs;
        }
      }
    }

    // data reordering:
    for(long m = 1, j = 0; m < my_n - 1; ++m)
    {
      for(long k = my_n >> 1; (!((j^=k)&k)); k>>=1)
      {
        ;
      }

      if(j > m)
      {
        std::swap(out[m], out[j]);
      }
    }

    // Normalize for backwards transform (done externally).
  }

  } // namespace detail

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_ALGORITHMS_HPP

