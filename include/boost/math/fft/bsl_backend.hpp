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
  #include <vector>
  #include <exception>
  #include <type_traits>

  #include <boost/math/fft/algorithms.hpp>
  #include <boost/math/fft/abstract_ring.hpp>

  namespace boost { namespace math {  namespace fft {


  /*
    Boost DFT backend:
    It handles RingTypes and it calls the appropriate specialized functions if
    the type is complex.
    
    Definition:
    A type is considered "complex" if it has ::value_type, and has sin() and cos()
    free functions.
    
    We have specialized algorithms for complex numbers and general purpose DFT
    that need the specification of a root of unity.
    The general purpose DFT work with complex and non-complex types, but its
    performance and precision could be lower than the specialized complex
    versions.
    
    This interface selects general purpose DFT for non-complex types,
    and for complex types the default behaviour is to use the specialized
    complex algorithms, unless the user provides a root of unity 'W', in which case
    the interface will execute general purpose DFT using W.
  */
  template<class RingType>
  class bsl_dft
  {
  private:
    enum plan_type { forward_plan , backward_plan};
    
    void execute(plan_type plan, const RingType * in, RingType* out)const
    {
      // select the implementation according to the type
      if(not has_root() and detail::is_complex<RingType>::value)
      {
          const int sign = (plan == forward_plan ? -1 : 1);
          // select the implementation according to the DFT size
          if( detail::is_power2(size())  )
            detail::dft_power2_dif(in,in+size(),out,sign);
          else
            detail::dft_complex_prime_bruteForce(in,in+size(),out,sign);
      }else
      {
        const RingType w = root();
        const RingType w_inverse = detail::power(w,size()-1);
        
        const RingType w_execute = (plan==forward_plan ? w : w_inverse);
        
        // select the implementation according to the DFT size
        if( detail::is_power2(size())  )
          detail::dft_power2_dit(in,in+size(),out,w_execute);
        else
          detail::dft_composite_dit(in,in+size(),out,w_execute);
      }
      
    }
    
  public:
    constexpr bsl_dft(std::size_t n=0)
      : my_size { n }
    { }
    
    // the provided root of unity is used instead of exp(-i 2 pi/n)
    constexpr bsl_dft(std::size_t n, RingType /* root of unity = */ w)
      : my_size { n }, _has_root{true}, my_root{w} , my_inverse_root{detail::power(w,n-1)}
    { }

    ~bsl_dft()
    {
    }
    
    bool has_root()const
    {
      return _has_root;
    }
    
    void resize(std::size_t new_size)
    {
      my_size = new_size;
      _has_root = false;
    }
    void resize(std::size_t new_size, RingType w)
    {
      my_size = new_size;
      _has_root = true;
      my_root = w;
      my_inverse_root = detail::power(w,new_size-1);
    }
    
    // non complex types
    template<typename U = RingType>
    RingType root(typename std::enable_if< detail::is_complex<U>::value == false >::type* = nullptr) const
    {
      if(not has_root())
        std::runtime_error("no root has been defined for this DFT size");
      return my_root;
    }
    // complex types
    template<typename U = RingType>
    RingType root(typename std::enable_if< detail::is_complex<U>::value == true >::type* = nullptr) const
    {
      if(has_root())
        return my_root;
      using std::sin;
      using std::cos;
      using real_value_type = typename RingType::value_type;
      real_value_type phase = 2*boost::math::constants::pi<real_value_type>()/size();
      return RingType(cos(phase),-sin(phase));
    }
    constexpr std::size_t size() const { return my_size; }

    void forward(const RingType* in, RingType* out) const
    {
      execute(forward_plan,in,out);   
    }

    void backward(const RingType* in, RingType* out) const
    {
      execute(backward_plan,in,out);   
    }

  private:
    std::size_t my_size{};
    bool _has_root=false;
    RingType my_root, my_inverse_root;
  };

  } } } // namespace boost::math::fft

#endif // BOOST_MATH_FFT_BSLBACKEND_HPP
