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


  namespace boost { namespace math {  namespace fft {

  namespace detail {

  // The const_helper (below) allows using constexpr whenever possible.
  // The allow_constexpr decides about whether constexpr is possible.
  // TODO: decide if we can find a better implementation for this.
  template <class T>
  constexpr bool allow_constexpr = (std::numeric_limits<T>::digits10 <= 15 and std::numeric_limits<T>::digits10 >= 6);

  template <class T> struct constexpr_consts {
    static constexpr T pi = boost::math::constants::pi<T>();
    static constexpr T _1 = T{1};
  };
  template<class T> constexpr T constexpr_consts<T>::pi;
  template<class T> constexpr T constexpr_consts<T>::_1;

  template <class T> struct const_consts {
    static const T pi;
    static const T _1;
  };
  template<class T> const T const_consts<T>::pi = boost::math::constants::pi<T>();
  template<class T> const T const_consts<T>::_1 = T{1};

  template <class T>
  using const_helper = typename std::conditional<allow_constexpr<T>, constexpr_consts<T>, const_consts<T>>::type;


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
        throw std::domain_error("power(x,n) expects n>1.");

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
  
  long least_power2(long x)
  {
    /*
      Returns the biggest power of two that is smaller or equal than x.
    */
    for(long lsb = x & -x; lsb!=x;)
    {
      x-=lsb;
      lsb = x&-x;
    }
    return x;
  }
  
  bool is_power2(long x) { return x == (x & -x);}
  
  template<class T>
  void dft_generic_prime_bruteForce(const T* in_first, const T* in_last, T* out, const T w)
  /*
    assumptions: 
    - allocated memory in out is enough to hold distance(in_first,in_last) element,
  */
  {
    const long N = std::distance(in_first,in_last);
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
  void dft_complex_prime_bruteForce(
    const complex_value_type* in_first, 
    const complex_value_type* in_last, 
    complex_value_type* out, int sign)
  /*
    assumptions: 
    - allocated memory in out is enough to hold distance(in_first,in_last) element,
  */
  {
    using real_value_type = typename complex_value_type::value_type;
    const long N = std::distance(in_first,in_last);
    const real_value_type inv_N = real_value_type{1}/N;
    const real_value_type signed_pi = sign * const_helper<real_value_type>::pi;
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
        real_value_type phase = 2*signed_pi * i *j * inv_N;

        using std::cos;
        using std::sin;

        sum += in_first[j] * complex_value_type{cos(phase),sin(phase)};
      }
      work_space[i] = sum;
    }
    
    if(out != work_space)
      std::copy(work_space,work_space+N,out);
  }
  
  template <class T>
  void dft_power2_dit(const T *in_first, const T *in_last, T* out, const T e)
  {
    /*
      Cooley-Tukey mapping, in-place Decimation in Time 
    */
    const long ptrdiff = std::distance(in_first,in_last);
    if(ptrdiff <=0 )
      return;
    const long n = least_power2(ptrdiff);
    
    if(in_first!=out)
      std::copy(in_first,in_last,out);
    
    if (n == 1)
        return;

    auto _1 = const_helper<T>::_1;
    
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
        T ej = _1;
        for (int j = 0; j < len / 2; ++j)
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
  void dft_power2_dif(
    const complex_value_type *in_first, 
    const complex_value_type *in_last, 
    complex_value_type* out, int sign)
  {
    // Naive in-place complex DFT.
    const long ptrdiff = std::distance(in_first,in_last);
    if(ptrdiff <=0 )
      return;
    const long my_n = least_power2(ptrdiff);
    
    if(in_first!=out)
      std::copy(in_first,in_last,out);
    
    if (my_n == 1)
        return;
    
    if(in_first!=out)
      std::copy(in_first, in_last, out);

    using local_float_type = typename complex_value_type::value_type;
    const local_float_type signed_pi = sign * const_helper<local_float_type>::pi;
    // Recursive decimation in frequency.
    for(long m = my_n; m > 1; m /= 2)
    {
      long mh = m / 2;

      const local_float_type phi = signed_pi / mh;

      for(long j = 0; j < mh; ++j)
      {
        const local_float_type p = phi * j;

        using std::cos;
        using std::sin;

        complex_value_type cs(cos(p), sin(p));

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

