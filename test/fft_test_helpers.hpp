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
class test_complex_dft_prime_rader
{
  /*
    Special backend for testing the complex_dft_prime_rader implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_complex_dft_prime_rader(std::size_t n)
    : my_size{n}
  { }
  
  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::complex_dft_prime_rader(in,in+size(),out,1);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::complex_dft_prime_rader(in,in+size(),out,-1);
  }
private:
  std::size_t my_size;
};

template<class NativeComplexType>
class test_dft_prime_bruteForce
{
  /*
    Special backend for testing the dft_prime_bruteForce implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_dft_prime_bruteForce(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    NativeComplexType w{detail::complex_root_of_unity<NativeComplexType>(size())};
    detail::dft_prime_bruteForce(in,in+size(),out,w);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    NativeComplexType w{detail::complex_inverse_root_of_unity<NativeComplexType>(size())};
    detail::dft_prime_bruteForce(in,in+size(),out,w);
  }

private:
  std::size_t my_size;
};

template<class NativeComplexType>
class test_complex_dft_prime_bruteForce
{
  /*
    Special backend for testing the complex_dft_prime_bruteForce implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_complex_dft_prime_bruteForce(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::complex_dft_prime_bruteForce(in,in+size(),out,1);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::complex_dft_prime_bruteForce(in,in+size(),out,-1);
  }

private:
  std::size_t my_size;
};

template<class NativeComplexType>
class test_dft_composite
{
  /*
    Special backend for testing the dft_composite
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_dft_composite(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    NativeComplexType w{detail::complex_root_of_unity<NativeComplexType>(size())};
    detail::dft_composite(in,in+size(),out,w);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    NativeComplexType w{detail::complex_inverse_root_of_unity<NativeComplexType>(size())};
    detail::dft_composite(in,in+size(),out,w);
  }

private:
  std::size_t my_size;
};

template<class NativeComplexType>
class test_complex_dft_composite
{
  /*
    Special backend for testing the complex_dft_composite
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_complex_dft_composite(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::complex_dft_composite(in,in+size(),out,1);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::complex_dft_composite(in,in+size(),out,-1);
  }

private:
  std::size_t my_size;
};

  
template<class NativeComplexType>
class test_dft_power2
{
  /*
    Special backend for testing the dft_power2 implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_dft_power2(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    NativeComplexType w{detail::complex_inverse_root_of_unity<NativeComplexType>(size())};
    detail::dft_power2(in,in+size(),out,w);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    NativeComplexType w{detail::complex_root_of_unity<NativeComplexType>(size())};
    detail::dft_power2(in,in+size(),out,w);
  }

private:
  std::size_t my_size;
};

template<class NativeComplexType>
class test_complex_dft_power2
{
  /*
    Special backend for testing the complex_dft_power2 implementation
  */
  using real_value_type = typename NativeComplexType::value_type;
public:
  constexpr test_complex_dft_power2(std::size_t n)
    : my_size{n}
  { }

  void resize(std::size_t new_size)
  {
    my_size = new_size;
  }
  constexpr std::size_t size() const { return my_size; }

  void forward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::complex_dft_power2(in,in+size(),out,1);
  }

  void backward(const NativeComplexType* in, NativeComplexType* out) const
  {
    detail::complex_dft_power2(in,in+size(),out,-1);
  }

private:
  std::size_t my_size;
};

namespace my_modulo_lib
{
  template <typename T, T x>
  class field_modulo
  {
   public:
    typedef T integer;
    static constexpr T mod = x;
  };

  /*
      Modular Integers
  */
  template <typename Field>
  class mint;

  template <typename Field>
  std::ostream& operator<<(std::ostream& os, const mint<Field>& A)
  {
    return os << A.x << " (mod " << Field::mod << ")";
  }
  template <typename Field>
  class mint
  {
    typedef typename Field::integer integer;
    integer x;

   public:
    constexpr mint() : x{0} {}

    template <typename int_type>
    mint(int_type _x) : x{_x}
    {
      x %= Field::mod;
      if (x < 0)
        x += Field::mod;
    }

    mint(const mint<Field>& that) : x{that.x} {}

    mint<Field>& operator=(const mint<Field>& that)
    {
      return x = that.x, *this;
    }

    // explicit operator bool() const { return x == integer{0}; }
    operator integer() const { return x; }

    mint<Field>& operator+=(const mint<Field>& that)
    {
      return x = (x + that.x) % Field::mod, *this;
    }

    mint<Field>& operator*=(const mint<Field>& t)
    {
      // direct multiplication
      x = (x * t.x) % Field::mod;
      return *this;
    }
    bool operator==(const mint<Field>& that) const { return x == that.x; }
    bool operator!=(const mint<Field>& that) const { return x != that.x; }

    mint<Field> inverse() const { return ::boost::math::fft::detail::power(*this, Field::mod - 2); }

    friend std::ostream& operator<<<Field>(std::ostream& os,
                                           const mint<Field>& A);
  };

  template <typename T>
  mint<T> operator+(const mint<T>& A, const mint<T>& B)
  {
    mint<T> C{A};
    return C += B;
  }
  template <typename T>
  mint<T> operator*(const mint<T>& A, const mint<T>& B)
  {
    mint<T> C{A};
    return C *= B;
  }
  template <typename T>
  mint<T>& operator/=(mint<T>& A, const mint<T>& B)
  {
    return A *= B.inverse();
  }
  template <typename T>
  mint<T> operator/(const mint<T>& A, const mint<T>& B)
  {
    mint<T> C{A};
    return C /= B;
  }
}  // namespace my_modulo_lib

}}} // namespace boost::math::fft
  
#endif // BOOST_MATH_FFT_TEST_HELPERS_HPP
