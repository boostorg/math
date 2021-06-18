///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FFT_BSLBACKEND_HPP
  #define BOOST_MATH_FFT_BSLBACKEND_HPP

  #include <algorithm>
  #include <cmath>
  #include <complex>
  #include <vector>

  #include <boost/math/constants/constants.hpp>

  namespace boost { namespace math {  namespace fft {

  namespace detail {
  
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
  template <class iter, class T>
  void generic_dft_power2(iter first, iter last, const T e)
  {
    const int n = std::distance(first, last);
    // TODO: determine power of two in a different way
    if (__builtin_popcount(n) != 1)
    // TODO: decide which exceptions to throw
        throw std::domain_error(std::string(__func__) +
                                 " n=" + std::to_string(n) +
                                 " must be a power of 2");
    if (n == 1)
        return;

    constexpr T _1 = T{1};
    
    int nbits = 0;
    std::vector<T> e2{e};
    for (int m = n / 2; m > 0; m >>= 1, ++nbits)
      e2.push_back(e2.back() * e2.back());

    std::reverse(e2.begin(), e2.end());

    iter _i = first;
    for (int i = 0; i < n; ++i, ++_i)
    {
      int ib = i, j = 0;
      iter _j = first;

      for (int b = 0; b < nbits; ib >>= 1, ++b)
          j = (j << 1) | (ib & 1);

      std::advance(_j, j);

      if (i < j)
          std::swap(*_i, *_j);
    }
    for (int len = 2, k = 1; len <= n; len <<= 1, ++k)
    {
      for (int i = 0; i < n; i += len)
      {
        T ej = _1;
        for (int j = 0; j < len / 2; ++j)
        {
          iter u = first + i + j, v = first + i + j + len / 2;
          T Bu = *u, Bv = *v * ej;
          *u = Bu + Bv;
          *v = Bu - Bv;
          ej *= e2[k];
        }
      }
    }
  }

  } // namespace detail

  template<class NativeComplexType>
  class bsl_dft
  {
  private:
    using real_value_type    = typename NativeComplexType::value_type;
    using complex_value_type = typename std::complex<real_value_type>;
    enum plan_type { forward_plan , backward_plan};
    
    void execute(plan_type plan, complex_value_type* out)const
    {
      // Naive in-place complex DFT.

      using local_float_type = typename complex_value_type::value_type;

      const local_float_type pi = (plan==backward_plan ?  boost::math::constants::pi<local_float_type>()
                                                       : -boost::math::constants::pi<local_float_type>());
      const long my_n = size();
      // Recursive decimation in frequency.
      for(long m = my_n; m > 1; m /= 2)
      {
        long mh = m / 2;

        const local_float_type phi = pi / mh;

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

  public:
    constexpr bsl_dft(std::ptrdiff_t n)
      : my_size          { n }
    { }

    ~bsl_dft()
    {
    }

    constexpr std::ptrdiff_t size() const { return my_size; }

    void forward(const complex_value_type* in, complex_value_type* out) const
    {
      if(in!=out)
        std::copy(in, in + size(), out);
      execute(forward_plan,out);   
    }

    void backward(const complex_value_type* in, complex_value_type* out) const
    {
      if(in!=out)
        std::copy(in, in + size(), out);
      execute(backward_plan,out);   
    }

  private:
    const std::ptrdiff_t      my_size;
  };
  
  template<class NativeComplexType>
  class generic_bsl_dft
  {
    using real_value_type = typename NativeComplexType::value_type;
    static constexpr real_value_type pi = boost::math::constants::pi<real_value_type>();
  public:
    constexpr generic_bsl_dft(std::ptrdiff_t n)
      : my_size{n}
    { }

    ~generic_bsl_dft()
    {
    }

    constexpr std::ptrdiff_t size() const { return my_size; }

    void forward(const NativeComplexType* in, NativeComplexType* out) const
    {
      if(in!=out)
        std::copy(in, in + size(), out);
      NativeComplexType w{std::cos(2*pi/size()),-std::sin(2*pi/size())};
      detail::generic_dft_power2(out,out+size(),w);
    }

    void backward(const NativeComplexType* in, NativeComplexType* out) const
    {
      if(in!=out)
        std::copy(in, in + size(), out);
      NativeComplexType w{std::cos(2*pi/size()),std::sin(2*pi/size())};
      detail::generic_dft_power2(out,out+size(),w);
    }

  private:
    const std::ptrdiff_t my_size;
  };

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_BSLBACKEND_HPP
