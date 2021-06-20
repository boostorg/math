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

  #include <boost/math/fft/algorithms.hpp>

  namespace boost { namespace math {  namespace fft {

  template<class NativeComplexType>
  class bsl_dft
  {
  private:
    using real_value_type    = typename NativeComplexType::value_type;
    using complex_value_type = typename std::complex<real_value_type>;
    enum plan_type { forward_plan , backward_plan};
    static constexpr real_value_type pi = boost::math::constants::pi<real_value_type>();
    
    void execute(plan_type plan, const complex_value_type * in, complex_value_type* out)const
    {
      const int sign = (plan == forward_plan ? -1 : 1);
      // select the implementation according to the DFT size
      if( detail::is_power2(size())  )
        detail::dft_power2_dif(in,in+size(),out,sign);
      else
      { 
        // NativeComplexType w{std::cos(2*pi/size()), sign * std::sin(2*pi/size())};
        detail::dft_complex_prime_bruteForce(in,in+size(),out,sign);
      }
    }
    
  public:
    constexpr bsl_dft(std::size_t n)
      : my_size { n }
    { }

    ~bsl_dft()
    {
    }
    
    void resize(std::size_t new_size)
    {
      my_size = new_size;
    }
    constexpr std::size_t size() const { return my_size; }

    void forward(const complex_value_type* in, complex_value_type* out) const
    {
      execute(forward_plan,in,out);   
    }

    void backward(const complex_value_type* in, complex_value_type* out) const
    {
      execute(backward_plan,in,out);   
    }

  private:
    std::size_t my_size;
  };

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_BSLBACKEND_HPP
