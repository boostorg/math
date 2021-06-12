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

  #include <boost/math/constants/constants.hpp>

  namespace boost { namespace math {  namespace fft {

  namespace detail {

  struct bsl_fft_plan
  {
    const int my_n;
    const int my_sign;

    constexpr bsl_fft_plan(const int n,
                           const int s = 1) : my_n   (n),
                                              my_sign(s) { }
  };

  template<typename T>
  struct bsl_traits
  {
    using plan_type = bsl_fft_plan;

    using real_value_type = T;

    using complex_value_type = std::complex<T>;

    static plan_type plan_construct(int n, int sign) { return bsl_fft_plan(n, sign); }

    static void plan_execute(plan_type plan, complex_value_type* out)
    {
      // Naive in-place complex DFT.

      using local_float_type = typename complex_value_type::value_type;

      const local_float_type pi = ((plan.my_sign >= 0) ?  boost::math::constants::pi<local_float_type>()
                                                       : -boost::math::constants::pi<local_float_type>());

      // Recursive decimation in frequency.
      for(long m = plan.my_n; m > 1; m /= 2)
      {
        long mh = m / 2;
        local_float_type phi = pi / mh;

        for(long j=0; j < mh; ++j)
        {
          const local_float_type p = phi * j;

          using std::cos;
          using std::sin;

          complex_value_type cs(cos(p), sin(p));

          for (long t1=j; t1 < j + plan.my_n; t1 += m)
          {
            complex_value_type u = out[t1];
            complex_value_type v = out[t1 + mh];

            out[t1]      =  u + v;
            out[t1 + mh] = (u - v) * cs;
          }
        }
      }

      // data reordering:
      for(long m = 1, j = 0; m < plan.my_n - 1; ++m)
      {
        for(long k = plan.my_n >> 1; (!((j^=k)&k)); k>>=1)
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

    static void plan_destroy(plan_type p) { ; }
  };

  } // namespace detail

  template<class NativeComplexType>
  class bsl_dft
  {
  private:
    using real_value_type    = typename NativeComplexType::value_type;
    using plan_type          = typename detail::bsl_traits<real_value_type>::plan_type;
    using complex_value_type = std::complex<real_value_type>;

  public:
    constexpr bsl_dft(std::ptrdiff_t n)
      : my_size          { n },
        my_forward_plan  { detail::bsl_traits<real_value_type>::plan_construct(size(), +1) },
        my_backward_plan { detail::bsl_traits<real_value_type>::plan_construct(size(), -1) }
    { }

    template<typename InputIteratorType,
             typename OutputIteratorType>
    bsl_dft(InputIteratorType  InFirst,
            InputIteratorType  InLast,
            OutputIteratorType OutFirst,
            typename std::enable_if<(   (std::is_convertible<InputIteratorType,  const complex_value_type*>::value == true)
                                     && (std::is_convertible<OutputIteratorType,       complex_value_type*>::value == true))>::type* = nullptr)
      : my_size          { std::distance(InFirst, InLast) },
        my_forward_plan  { detail::bsl_traits<real_value_type>::plan_construct(size(), +1) },
        my_backward_plan { detail::bsl_traits<real_value_type>::plan_construct(size(), -1) }
    { }

    ~bsl_dft()
    {
      detail::bsl_traits<real_value_type>::plan_destroy(my_forward_plan);
      detail::bsl_traits<real_value_type>::plan_destroy(my_backward_plan);
    }

    constexpr std::ptrdiff_t size() const { return my_size; }

    void forward(const complex_value_type* in, complex_value_type* out) const
    {
      using local_complex_type = typename detail::bsl_traits<real_value_type>::complex_value_type;

      std::copy(in, in + size(), out);

      detail::bsl_traits<real_value_type>::plan_execute
      (
        my_forward_plan,
        reinterpret_cast<local_complex_type*>(out)
      );
    }

    void backward(const complex_value_type* in, complex_value_type* out) const
    {
      using local_complex_type = typename detail::bsl_traits<real_value_type>::complex_value_type;

      std::copy(in, in + size(), out);

      detail::bsl_traits<real_value_type>::plan_execute
      (
        my_backward_plan,
        reinterpret_cast<local_complex_type*>(out)
      );
    }

  private:
    const std::ptrdiff_t my_size;
          plan_type      my_forward_plan;
          plan_type      my_backward_plan;
  };

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_BSLBACKEND_HPP
