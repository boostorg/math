///////////////////////////////////////////////////////////////////
//  Copyright Eduardo Quintana 2021
//  Copyright Janek Kozicki 2021
//  Copyright Christopher Kormanyos 2021
//  Distributed under the Boost Software License,
//  Version 1.0. (See accompanying file LICENSE_1_0.txt
//  or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_FFT_TEST_HELPERS_HPP
  #define BOOST_MATH_FFT_TEST_HELPERS_HPP

#include <boost/math/fft/algorithms.hpp>
  
namespace boost { namespace math { namespace fft {

template<class NativeComplexType>
class test_dft_generic_prime_bruteForce
{
  /*
    Special backend for testing the dft_prime_bruteForce implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_dft_generic_prime_bruteForce(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    using std::cos;
    using std::sin;
    using boost::math::constants::pi;

    NativeComplexType w{cos(2*pi<real_value_type>()/size()),-sin(2*pi<real_value_type>()/size())};
    detail::dft_generic_prime_bruteForce(in,in+size(),out,w);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    using std::cos;
    using std::sin;
    using boost::math::constants::pi;

    NativeComplexType w{cos(2*pi<real_value_type>()/size()),sin(2*pi<real_value_type>()/size())};
    detail::dft_generic_prime_bruteForce(in,in+size(),out,w);
  }

private:
  std::size_t my_size;
};
template<class NativeComplexType>
class test_dft_generic_composite
{
  /*
    Special backend for testing the dft_composite_dit
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_dft_generic_composite(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    using std::cos;
    using std::sin;
    using boost::math::constants::pi;

    NativeComplexType w{cos(2*pi<real_value_type>()/size()),-sin(2*pi<real_value_type>()/size())};
    detail::dft_composite_dit(in,in+size(),out,w);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    using std::cos;
    using std::sin;
    using boost::math::constants::pi;

    NativeComplexType w{cos(2*pi<real_value_type>()/size()),sin(2*pi<real_value_type>()/size())};
    detail::dft_composite_dit(in,in+size(),out,w);
  }

private:
  std::size_t my_size;
};

template<class NativeComplexType>
class test_dft_complex_prime_bruteForce
{
  /*
    Special backend for testing the dft_prime_bruteForce implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_dft_complex_prime_bruteForce(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::dft_complex_prime_bruteForce(in,in+size(),out,-1);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::dft_complex_prime_bruteForce(in,in+size(),out,1);
  }

private:
  std::size_t my_size;
};
  
template<class NativeComplexType>
class test_dft_power2_dit
{
  /*
    Special backend for testing the dft_power2_dit implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_dft_power2_dit(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    using std::cos;
    using std::sin;
    using boost::math::constants::pi;

    NativeComplexType w{cos(2*pi<real_value_type>()/size()),-sin(2*pi<real_value_type>()/size())};
    detail::dft_power2_dit(in,in+size(),out,w);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    using std::cos;
    using std::sin;
    using boost::math::constants::pi;

    NativeComplexType w{cos(2*pi<real_value_type>()/size()),sin(2*pi<real_value_type>()/size())};
    detail::dft_power2_dit(in,in+size(),out,w);
  }

private:
  std::size_t my_size;
};

template<class NativeComplexType>
class test_dft_power2_dif
{
  /*
    Special backend for testing the dft_power2_dif implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_dft_power2_dif(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::dft_power2_dif(in,in+size(),out,-1);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::dft_power2_dif(in,in+size(),out,1);
  }

private:
  std::size_t my_size;
};

}}} // namespace boost::math::fft
  
#endif // BOOST_MATH_FFT_TEST_HELPERS_HPP
