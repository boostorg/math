//  (C) Copyright Nick Thompson 2019.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_TRIGONOMETRIC_HPP
#define BOOST_MATH_INTERPOLATORS_DETAIL_CARDINAL_TRIGONOMETRIC_HPP
#ifndef BOOST_MATH_BUILD_MODULE
#include <cstddef>
#include <cmath>
#include <stdexcept>
#include <vector>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/fft.hpp>
#include <boost/math/tools/detail/complex_from_real.hpp>

BOOST_MATH_NAMESPACE_BEGIN namespace interpolators { namespace detail {

template<typename Real>
class cardinal_trigonometric_detail {
public:
  using complex_type = typename BOOST_MATH_NAMESPACE::tools::detail::complex_from_real<Real>::type;

  template<class RandomAccessContainer>
  cardinal_trigonometric_detail(RandomAccessContainer const & v, Real t0, Real h) : m_t0{t0}, m_h{h}
  {
    std::size_t length = v.size();
    if (length == 0)
    {
      throw std::logic_error("At least one sample is required.");
    }
    if (h <= 0)
    {
      throw std::logic_error("The step size must be > 0");
    }
    // The period must be stored, since it cannot be recovered from the number of Fourier coefficients:
    m_T = m_h*static_cast<Real>(length);
    m_gamma.resize(length/2 + 1);
    BOOST_MATH_NAMESPACE::tools::real_fft<complex_type> fft(length);
    fft.transform(v, m_gamma);

    Real denom = static_cast<Real>(length);
    for (auto& g : m_gamma)
    {
      g = complex_type(Real(g.real()/denom), Real(g.imag()/denom));
    }
    // The transform of real data has a real mean; drop any rounding in its imaginary part.
    m_gamma[0] = complex_type(Real(m_gamma[0].real()), Real(0));

    if (length % 2 == 0)
    {
      m_gamma.back() = complex_type(Real(m_gamma.back().real()/2), Real(m_gamma.back().imag()));
    }
  }

  Real operator()(Real t) const
  {
    using std::sin;
    using std::cos;
    using BOOST_MATH_NAMESPACE::constants::two_pi;
    Real s = m_gamma[0].real();
    Real x = two_pi<Real>()*(t - m_t0)/m_T;
    Real z0 = cos(x);
    Real z1 = sin(x);
    // Horner's method on the complex polynomial with b <- gamma[k] + b*z:
    Real b0 = 0;
    Real b1 = 0;
    for (std::size_t k = m_gamma.size() - 1; k >= 1; --k)
    {
      Real u0 = b0*z0 - b1*z1;
      Real u1 = b0*z1 + b1*z0;
      b0 = m_gamma[k].real() + u0;
      b1 = m_gamma[k].imag() + u1;
    }

    s += 2*(b0*z0 - b1*z1);
    return s;
  }

  Real prime(Real t) const
  {
    using std::sin;
    using std::cos;
    using BOOST_MATH_NAMESPACE::constants::two_pi;
    Real x = two_pi<Real>()*(t - m_t0)/m_T;
    Real z0 = cos(x);
    Real z1 = sin(x);
    Real b0 = 0;
    Real b1 = 0;
    for (std::size_t k = m_gamma.size() - 1; k >= 1; --k)
    {
      Real u0 = b0*z0 - b1*z1;
      Real u1 = b0*z1 + b1*z0;
      b0 = static_cast<Real>(k)*m_gamma[k].real() + u0;
      b1 = static_cast<Real>(k)*m_gamma[k].imag() + u1;
    }
    // b*z = (b0*z0 - b1*z1) + i(b1*z0 + b0*z1)
    return -2*two_pi<Real>()*(b1*z0 + b0*z1)/m_T;
  }

  Real double_prime(Real t) const
  {
    using std::sin;
    using std::cos;
    using BOOST_MATH_NAMESPACE::constants::two_pi;
    Real x = two_pi<Real>()*(t - m_t0)/m_T;
    Real z0 = cos(x);
    Real z1 = sin(x);
    Real b0 = 0;
    Real b1 = 0;
    for (std::size_t k = m_gamma.size() - 1; k >= 1; --k)
    {
      Real u0 = b0*z0 - b1*z1;
      Real u1 = b0*z1 + b1*z0;
      b0 = static_cast<Real>(k*k)*m_gamma[k].real() + u0;
      b1 = static_cast<Real>(k*k)*m_gamma[k].imag() + u1;
    }
    return -2*two_pi<Real>()*two_pi<Real>()*(b0*z0 - b1*z1)/(m_T*m_T);
  }

  Real period() const
  {
    return m_T;
  }

  Real integrate() const
  {
    return m_T*m_gamma[0].real();
  }

  Real squared_l2() const
  {
    Real s = 0;
    // Add the smallest terms first for accuracy.
    for (std::size_t i = m_gamma.size() - 1; i >= 1; --i)
    {
      s += (m_gamma[i].real()*m_gamma[i].real() + m_gamma[i].imag()*m_gamma[i].imag());
    }
    s *= 2;
    s += m_gamma[0].real()*m_gamma[0].real();
    return s*m_T;
  }

private:
  Real m_t0;
  Real m_h;
  Real m_T;
  std::vector<complex_type> m_gamma;
};

}} BOOST_MATH_NAMESPACE_END
#endif
