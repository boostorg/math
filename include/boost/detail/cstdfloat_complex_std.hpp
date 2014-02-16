///////////////////////////////////////////////////////////////////////////////
// Copyright Christopher Kormanyos 2014.
// Copyright John Maddock 2014.
// Copyright Paul Bristow 2014.
// Distributed under the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)
//

// Implement a specialization of std::complex<> for *anything* that
// is defined as BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE.

#ifndef _BOOST_CSTDFLOAT_COMPLEX_STD_2014_02_15_HPP_
  #define _BOOST_CSTDFLOAT_COMPLEX_STD_2014_02_15_HPP_

  #include <complex>
  #include <boost/math/constants/constants.hpp>

  namespace std
  {
    template<>
    class complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>
    {
    public:
      typedef BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE value_type;

      #if defined(BOOST_NO_CXX11_CONSTEXPR)

      complex(const complex<value_type>& z) : re(z.real()),
                                              im(z.imag()) { }

      complex(value_type r = value_type(),
              value_type i = value_type()) : re(r),
                                             im(i) { }

      explicit complex(const complex<float>&);
      explicit complex(const complex<double>&);
      explicit complex(const complex<long double>&);

      const value_type& real() const { return re; }
      const value_type& imag() const { return im; }
            value_type& real()       { return re; }
            value_type& imag()       { return im; }
      #else

      BOOST_CONSTEXPR complex(const complex<value_type>& z) : re(z.real()),
                                                              im(z.imag()) { }

      BOOST_CONSTEXPR complex(value_type r = value_type(),
                              value_type i = value_type()) : re(r),
                                                              im(i) { }

      explicit BOOST_CONSTEXPR complex(const complex<float>&);
      explicit BOOST_CONSTEXPR complex(const complex<double>&);
      explicit BOOST_CONSTEXPR complex(const complex<long double>&);

      BOOST_CONSTEXPR value_type real() { return re; }
      BOOST_CONSTEXPR value_type imag() { return im; }
      #endif

      void real(value_type v) { re = v; }
      void imag(value_type v) { im = v; }

      complex<value_type>& operator+=(value_type v)
      {
        re += v;
        return *this;
      }

      complex<value_type>& operator-=(value_type v)
      {
        re -= v;
        return *this;
      }

      complex<value_type>& operator*=(value_type v)
      {
        re *= v;
        im *= v;

        return *this;
      }

      complex<value_type>& operator/=(value_type v)
      {
        re /= v;
        im /= v;

        return *this;
      }

      template<typename X>
      complex<value_type>& operator+=(const complex<X>& x)
      {
        re += value_type(x.re);
        im += value_type(x.im);

        return *this;
      }

      template<class X>
      complex<value_type>& operator-=(const complex<X>& x)
      {
        re -= value_type(x.re);
        im -= value_type(x.im);

        return *this;
      }

      template<class X>
      complex<value_type>& operator*=(const complex<X>& x)
      {
        const value_type re_x(x.re);
        const value_type im_x(x.im);

        const value_type tmp_re((re * re_x) - (im * im_x));
        const value_type tmp_im((re * im_x) + (im * re_x));

        re = tmp_re;
        im = tmp_im;

        return *this;
      }

      template<class X>
      complex<value_type>& operator/=(const complex<X>& x)
      {
        using std::sqrt;

        const value_type re_x(x.re);
        const value_type im_x(x.im);

        const value_type the_norm = sqrt((re_x * re_x) + (im_x * im_x));

        const value_type tmp_re = ((re * re_x) + (im * im_x)) / the_norm;
        const value_type tmp_im = ((im * re_x) - (re * im_x)) / the_norm;

        re = tmp_re;
        im = tmp_im;

        return *this;
      }

      template<class X>
      complex<value_type>& operator=(const complex<X>& x)
      {
        re = x.real();
        im = x.imag();

        return *this;
      }

      complex<value_type>& operator=(const value_type& x)
      {
        re = x;
        im = value_type(0);

        return *this;
      }

    private:
      value_type re;
      value_type im;
    };

    // Constructors from built-in floating-point types.
    #if defined(BOOST_NO_CXX11_CONSTEXPR)
    complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::complex(const complex<float>&        f) : re(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE( f.real())), im(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE( f.imag())) { }
    complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::complex(const complex<double>&       d) : re(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE( d.real())), im(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE( d.imag())) { }
    complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::complex(const complex<long double>& ld) : re(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(ld.real())), im(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(ld.imag())) { }
    #else
    BOOST_CONSTEXPR complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::complex(const complex<float>&        f) : re(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE( f.real())), im(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE( f.imag())) { }
    BOOST_CONSTEXPR complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::complex(const complex<double>&       d) : re(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE( d.real())), im(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE( d.imag())) { }
    BOOST_CONSTEXPR complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::complex(const complex<long double>& ld) : re(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(ld.real())), im(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(ld.imag())) { }
    #endif

    // Forward template function declarations.
    #if defined(BOOST_NO_CXX11_CONSTEXPR)
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& real(complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& imag(complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    #else
    inline BOOST_CONSTEXPR BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE real(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline BOOST_CONSTEXPR BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE imag(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    #endif
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE abs (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE arg (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE norm(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> conj (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    #if !defined(BOOST_NO_CXX11_FULL_COMPLEX_TRANSCENDENTALS)
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> proj (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    #endif
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> polar(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE&,
                                                                      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE&);

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sqrt (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sin  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> cos  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> tan  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    #if !defined(BOOST_NO_CXX11_FULL_COMPLEX_TRANSCENDENTALS)
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> asin (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> acos (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> atan (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    #endif
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> exp  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> log  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> log10(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow  (const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE&, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sinh (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> cosh (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> tanh (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    #if !defined(BOOST_NO_CXX11_FULL_COMPLEX_TRANSCENDENTALS)
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> asinh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> acosh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> atanh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
    #endif

    template<class char_type, class traits_type>
    inline std::basic_ostream<char_type, traits_type>& operator<<(std::basic_ostream<char_type, traits_type>&, const std::complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);

    template<class char_type, class traits_type>
    inline std::basic_istream<char_type, traits_type>& operator>>(std::basic_istream<char_type, traits_type>&, std::complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>&);
  } // namespace std

  // Here is a simple helper function.
  namespace boost { namespace cstdfloat { namespace detail {
  template<class float_type> std::complex<float_type> iz_helper__x(const std::complex<float_type>& x)
  {
    const float_type tmp_r = x.real();
    return std::complex<float_type>(-x.imag(), tmp_r);
  }
  } } } // boost::cstdfloat::detail

  namespace std
  {
    // ISO/IEC 14882:2011, Section 26.4.7, specific values.
    #if defined(BOOST_NO_CXX11_CONSTEXPR)
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& real(complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x) { return x.real(); }
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& imag(complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x) { return x.imag(); }
    #else
    inline BOOST_CONSTEXPR BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE real(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x) { return x.real(); }
    inline BOOST_CONSTEXPR BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE imag(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x) { return x.imag(); }
    #endif
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE abs (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x) { using std::sqrt;  return sqrt ((real(x) * real(x)) + (imag(x) * imag(x))); }
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE arg (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x) { using std::atan2; return atan2(x.imag(), x.real()); }
    inline BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE norm(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x) { return (real(x) * real(x)) + (imag(x) * imag(x)); }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> conj (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(x.real(), -x.imag()); }
    #if !defined(BOOST_NO_CXX11_FULL_COMPLEX_TRANSCENDENTALS)
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> proj (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE two_over_denom = BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(2) / (std::norm(x) + 1);

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(x.real() * two_over_denom, x.imag() * two_over_denom);
    }
    #endif
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> polar(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& rho,
                                                                      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& theta)
    {
      using std::sin;
      using std::cos;

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(rho * cos(theta), rho * sin(theta));
    }

    // Global add, sub, mul, div.
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator+(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u.real() + v.real(), u.imag() + v.imag()); }
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator-(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u.real() - v.real(), u.imag() - v.imag()); }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator*(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v)
    {
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::value_type ur = u.real();
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::value_type ui = u.imag();
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::value_type vr = v.real();
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>::value_type vi = v.imag();

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>((ur * vr) - (ui * vi), (ur * vi) + (ui * vr));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator/(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v)
    {
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE the_norm = std::norm(v);

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(((u.real() * v.real()) + (u.imag() * v.imag())) / the_norm,
                                                                  ((u.imag() * v.real()) - (u.real() * v.imag())) / the_norm);
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator+(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u.real() + v, u.imag()); }
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator-(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u.real() - v, u.imag()); }
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator*(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u.real() * v, u.imag() * v); }
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator/(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u.real() / v, u.imag() / v); }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator+(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u + v.real(), v.imag()); }
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator-(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u - v.real(), -v.imag()); }
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator*(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(u * v.real(), u * v.imag()); }
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator/(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE v_norm = norm(v); return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>((u * v.real()) / v_norm, (-u * v.imag()) / v_norm); }

    // Unary plus / minus.
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator+(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u) { return u; }
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> operator-(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u) { return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(-u.real(), -u.imag()); }

    // Equality and inequality.
    #if defined(BOOST_NO_CXX11_CONSTEXPR)
    inline bool operator==(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return ((u.real() == v.real()) && (u.imag() == v.imag())); }
    inline bool operator!=(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return ((u.real() != v.real()) || (u.imag() != v.imag())); }

    inline bool operator==(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& v) { return ((u.real() == v) && (u.imag() == BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(0))); }
    inline bool operator!=(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& v) { return ((u.real() != v) || (u.imag() != BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(0))); }

    inline bool operator==(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return ((u == v.real()) && (v.imag() == BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(0))); }
    inline bool operator!=(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return ((u != v.real()) || (v.imag() != BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(0))); }
    #else
    inline BOOST_CONSTEXPR bool operator==(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return ((u.real() == v.real()) && (u.imag() == v.imag())); }
    inline BOOST_CONSTEXPR bool operator!=(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return ((u.real() != v.real()) || (u.imag() != v.imag())); }

    inline BOOST_CONSTEXPR bool operator==(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& v) { return ((u.real() == v) && (u.imag() == BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(0))); }
    inline BOOST_CONSTEXPR bool operator!=(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& u, const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& v) { return ((u.real() != v) || (u.imag() != BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(0))); }

    inline BOOST_CONSTEXPR bool operator==(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return ((u == v.real()) && (v.imag() == BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(0))); }
    inline BOOST_CONSTEXPR bool operator!=(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& u, const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& v) { return ((u != v.real()) || (v.imag() != BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(0))); }
    #endif

    // ISO/IEC 14882:2011, Section 26.4.8, transcendentals.
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sqrt (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sin  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> cos  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> tan  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> asin (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> acos (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> atan (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> exp  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> log  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> log10(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x,
                                                                      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& a);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow  (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x,
                                                                      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& a);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow  (const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& x,
                                                                      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& a);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sinh (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> cosh (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> tanh (const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> asinh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> acosh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);
    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> atanh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x);

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sqrt(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      using std::fabs;
      using std::sqrt;

      // sqrt(x) for (x in C) =
      //  (s       , xi / 2s) : for xr > 0,
      //  (|xi| / 2s, +-s)    : for xr < 0,
      //  (sqrt(xi), sqrt(xi) : for xr = 0,
      // where s = sqrt{ [ |xr| + sqrt(xr^2 + xi^2) ] / 2 },
      // and the +- sign is the same as the sign of xi.

      if(x.real() > 0)
      {
        const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE s = sqrt((fabs(x.real()) + std::abs(x)) / 2);

        return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(s, x.imag() / (s * 2));
      }
      else if(x.real() < 0)
      {
        const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE s = sqrt((fabs(x.real()) + std::abs(x)) / 2);

        const bool imag_is_neg = (x.imag() < 0);

        return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(fabs(x.imag()) / (s * 2), (imag_is_neg ? -s : s));
      }
      else
      {
        const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sqrt_xi_half = sqrt(x.imag() / 2);

        return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(sqrt_xi_half, sqrt_xi_half);
      }
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sin(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      using std::sin;
      using std::cos;
      using std::exp;

      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sin_x  = sin (x.real());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cos_x  = cos (x.real());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_yp = exp (x.imag());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_ym = BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) / exp_yp;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sinh_y = (exp_yp - exp_ym) / 2;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cosh_y = (exp_yp + exp_ym) / 2;

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(sin_x * cosh_y, cos_x * sinh_y);
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> cos(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      using std::sin;
      using std::cos;
      using std::exp;

      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sin_x  = sin (x.real());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cos_x  = cos (x.real());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_yp = exp (x.imag());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_ym = BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) / exp_yp;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sinh_y = (exp_yp - exp_ym) / 2;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cosh_y = (exp_yp + exp_ym) / 2;

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(cos_x * cosh_y, -(sin_x * sinh_y));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> tan(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      using std::sin;
      using std::cos;
      using std::exp;

      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sin_x  = sin (x.real());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cos_x  = cos (x.real());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_yp = exp (x.imag());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_ym = BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) / exp_yp;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sinh_y = (exp_yp - exp_ym) / 2;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cosh_y = (exp_yp + exp_ym) / 2;

      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE a =   sin_x * cosh_y;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE b =   cos_x * sinh_y;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE c =   cos_x * cosh_y;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE d = -(sin_x * sinh_y);

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>((a * c) + (b * d), (b * c) - (a * d)) / ((c * c) + (d * d));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> asin(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      return -boost::cstdfloat::detail::iz_helper__x(std::log(boost::cstdfloat::detail::iz_helper__x(x) + std::sqrt(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) - (x * x))));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> acos(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      return boost::math::constants::half_pi<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>() - std::asin(x);
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> atan(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> izz = boost::cstdfloat::detail::iz_helper__x(x);

      return boost::cstdfloat::detail::iz_helper__x(std::log(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) - izz) - std::log(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) + izz)) / 2;
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> exp(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      using std::exp;

      return std::polar(exp(x.real()), x.imag());
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> log(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      using std::atan2;
      using std::log;

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(log(std::norm(x)) / 2, atan2(x.imag(), x.real()));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> log10(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      return std::log(x) / boost::math::constants::ln_ten<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>();
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x,
                                                                    const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& a)
    {
      std::exp(a * std::log(x));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x,
                                                                    const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& a)
    {
      std::exp(a * std::log(x));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> pow(const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE& x,
                                                                    const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& a)
    {
      std::exp(a * std::log(x));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> sinh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      using std::sin;
      using std::cos;
      using std::exp;

      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sin_y  = sin (x.imag());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cos_y  = cos (x.imag());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_xp = exp (x.real());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_xm = BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) / exp_xp;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sinh_x = (exp_xp - exp_xm) / 2;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cosh_x = (exp_xp + exp_xm) / 2;

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(cos_y * sinh_x, cosh_x * sin_y);
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> cosh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      using std::sin;
      using std::cos;
      using std::exp;

      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sin_y  = sin (x.imag());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cos_y  = cos (x.imag());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_xp = exp (x.real());
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE exp_xm = BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) / exp_xp;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE sinh_x = (exp_xp - exp_xm) / 2;
      const BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE cosh_x = (exp_xp + exp_xm) / 2;

      return complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(cos_y * cosh_x, sin_y * sinh_x);
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> tanh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> ex_plus  = std::exp(x);
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> ex_minus = 1 / ex_plus;

      return (ex_plus - ex_minus) / (ex_plus + ex_minus);
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> asinh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      return std::log(x + std::sqrt((x * x) + BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1)));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> acosh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> zp(x.real() + 1, x.imag());
      const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> zm(x.real() - 1, x.imag());

      return std::log(x + (zp * std::sqrt(zm / zp)));
    }

    inline complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE> atanh(const complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      return (std::log(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) + x) - std::log(BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE(1) - x)) / 2.0;
    }

    template<class char_type, class traits_type>
    inline std::basic_ostream<char_type, traits_type>& operator<<(std::basic_ostream<char_type, traits_type>& os, const std::complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      std::basic_ostringstream<char_type, traits_type> ostr;
      ostr.flags(os.flags());
      ostr.imbue(os.getloc());
      ostr.precision(os.precision());

      ostr << char_type('(')
            << x.real()
            << char_type(',')
            << x.imag()
            << char_type(')');

      return (os << ostr.str());
    }

    template<class char_type, class traits_type>
    inline std::basic_istream<char_type, traits_type>& operator>>(std::basic_istream<char_type, traits_type>& is, std::complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>& x)
    {
      BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE rx;
      BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE ix;

      char_type the_char;

      static_cast<void>(is >> the_char);

      if(the_char == static_cast<char_type>('('))
      {
        static_cast<void>(is >> rx >> the_char);

        if(the_char == static_cast<char_type>(','))
        {
          static_cast<void>(is >> ix >> the_char);

          if(the_char == static_cast<char_type>(')'))
          {
            x = complex<BOOST_CSTDFLOAT_EXTENDED_COMPLEX_FLOAT_TYPE>(rx, ix);
          }
          else
          {
            is.setstate(ios_base::failbit);
          }
        }
        else if(the_char == static_cast<char_type>(')'))
        {
          x = rx;
        }
        else
        {
          is.setstate(ios_base::failbit);
        }
      }
      else
      {
        static_cast<void>(is.putback(the_char));

        static_cast<void>(is >> rx);

        x = rx;
      }

      return is;
    }
  } // namespace std

#endif // _BOOST_CSTDFLOAT_COMPLEX_STD_2014_02_15_HPP_
