//  (C) Copyright Nick Thompson 2017.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_CHEBYSHEV_TRANSFORM_HPP
#define BOOST_MATH_SPECIAL_CHEBYSHEV_TRANSFORM_HPP
#ifndef BOOST_MATH_BUILD_MODULE
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>
#endif
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/chebyshev.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/tools/fft.hpp>
#include <boost/math/tools/detail/complex_from_real.hpp>

BOOST_MATH_NAMESPACE_BEGIN

namespace detail{

// The type II discrete cosine transform with FFTW's REDFT10 normalisation,
//
//    y[k] = 2 sum_{j=0}^{n-1} x[j] cos(pi*k*(2j+1)/(2n)),
//
// computed with one real FFT of length n (Makhoul, IEEE Trans. ASSP 28(1), 1980):
// with v[j] = x[2j] and v[n-1-j] = x[2j+1], and V the DFT of v,
// y[k] = 2 Re(exp(-i*pi*k/(2n)) V[k]) and y[n-k] = -2 Im(exp(-i*pi*k/(2n)) V[k]).
template <class Real>
class dct_ii
{
public:
   using complex_type = typename BOOST_MATH_NAMESPACE::tools::detail::complex_from_real<Real>::type;

   explicit dct_ii(std::size_t n) : m_n(n), m_fft(n), m_v(n), m_V(n / 2 + 1), m_cos(n / 2 + 1), m_sin(n / 2 + 1)
   {
      // exp(-i*pi*k/(2n)) = cos - i sin for k in [0, n/2]; the argument k/(2n) has at most one rounding.
      for (std::size_t k = 0; k <= n / 2; ++k)
      {
         Real x = static_cast<Real>(k) / static_cast<Real>(2 * n);
         m_cos[k] = BOOST_MATH_NAMESPACE::cos_pi(x);
         m_sin[k] = BOOST_MATH_NAMESPACE::sin_pi(x);
      }
   }

   // In place: replaces the n values at x with their transform.
   void execute(Real* x)
   {
      execute(x, x);
   }

   // Reads n values from x and writes n values to y. x and y are either the same array (all of x is read
   // into the workspace before y is written) or don't overlap. Any n > 0 works, odd or even; the constructor
   // allocates all the workspace, so execute never allocates.
   void execute(const Real* x, Real* y)
   {
      for (std::size_t j = 0; 2 * j < m_n; ++j)
      {
         m_v[j] = x[2 * j];
      }
      for (std::size_t j = 0; 2 * j + 1 < m_n; ++j)
      {
         m_v[m_n - 1 - j] = x[2 * j + 1];
      }
      m_fft.transform(m_v, m_V);
      for (std::size_t k = 0; k <= m_n / 2; ++k)
      {
         Real vr = m_V[k].real();
         Real vi = m_V[k].imag();
         // (cos - i sin)(vr + i vi)
         y[k] = 2 * (m_cos[k] * vr + m_sin[k] * vi);
         if (k != 0 && m_n - k != k)
         {
            y[m_n - k] = -2 * (m_cos[k] * vi - m_sin[k] * vr);
         }
      }
   }

private:
   std::size_t m_n;
   BOOST_MATH_NAMESPACE::tools::real_fft<complex_type> m_fft;
   std::vector<Real> m_v;
   std::vector<complex_type> m_V;
   std::vector<Real> m_cos;
   std::vector<Real> m_sin;
};

}

BOOST_MATH_EXPORT template<class Real>
class chebyshev_transform
{
public:
    template<class F>
    chebyshev_transform(const F& f, Real a, Real b,
       Real tol = 500 * std::numeric_limits<Real>::epsilon(),
       size_t max_refinements = 16) : m_a(a), m_b(b)
    {
        if (a >= b)
        {
            throw std::domain_error("a < b is required.\n");
        }
        using BOOST_MATH_NAMESPACE::constants::half;
        using std::abs;
        Real bma = (b-a)*half<Real>();
        Real bpa = (b+a)*half<Real>();
        size_t n = 256;
        std::vector<Real> vf;

        size_t refinements = 0;
        while(refinements < max_refinements)
        {
            vf.resize(n);
            m_coeffs.resize(n);

            detail::dct_ii<Real> dct(n);
            Real inv_n = 1/static_cast<Real>(n);
            for(size_t j = 0; j < n/2; ++j)
            {
                // Use symmetry cos((j+1/2)pi/n) = - cos((n-1-j+1/2)pi/n)
                Real y = BOOST_MATH_NAMESPACE::cos_pi(static_cast<Real>(2*j+1)/static_cast<Real>(2*n));
                vf[j] = f(y*bma + bpa)*inv_n;
                vf[n-1-j]= f(bpa-y*bma)*inv_n;
            }

            dct.execute(vf.data(), m_coeffs.data());
            Real max_coeff = 0;
            for (auto const & coeff : m_coeffs)
            {
                if (abs(coeff) > max_coeff)
                {
                    max_coeff = abs(coeff);
                }
            }
            size_t j = m_coeffs.size() - 1;
            while (abs(m_coeffs[j])/max_coeff < tol)
            {
                --j;
            }
            // If ten coefficients are eliminated, the we say we've done all
            // we need to do:
            if (n - j > 10)
            {
                m_coeffs.resize(j+1);
                return;
            }

            n *= 2;
            ++refinements;
        }
    }

    inline Real operator()(Real x) const
    {
        return chebyshev_clenshaw_recurrence(m_coeffs.data(), m_coeffs.size(), m_a, m_b, x);
    }

    // Integral over entire domain [a, b]
    Real integrate() const
    {
          Real Q = m_coeffs[0]/2;
          for(size_t j = 2; j < m_coeffs.size(); j += 2)
          {
              Q += -m_coeffs[j]/((j+1)*(j-1));
          }
          return (m_b - m_a)*Q;
    }

    const std::vector<Real>& coefficients() const
    {
        return m_coeffs;
    }

    Real prime(Real x) const
    {
        Real z = (2*x - m_a - m_b)/(m_b - m_a);
        Real dzdx = 2/(m_b - m_a);
        if (m_coeffs.size() < 2)
        {
            return 0;
        }
        Real b2 = 0;
        Real d2 = 0;
        Real b1 = m_coeffs[m_coeffs.size() -1];
        Real d1 = 0;
        for(size_t j = m_coeffs.size() - 2; j >= 1; --j)
        {
            Real tmp1 = 2*z*b1 - b2 + m_coeffs[j];
            Real tmp2 = 2*z*d1 - d2 + 2*b1;
            b2 = b1;
            b1 = tmp1;

            d2 = d1;
            d1 = tmp2;
        }
        return dzdx*(z*d1 - d2 + b1);
    }

private:
    std::vector<Real> m_coeffs;
    Real m_a;
    Real m_b;
};

BOOST_MATH_NAMESPACE_END
#endif
