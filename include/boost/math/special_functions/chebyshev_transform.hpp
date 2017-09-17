//  (C) Copyright Nick Thompson 2017.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_CHEBYSHEV_TRANSFORM_HPP
#define BOOST_MATH_SPECIAL_CHEBYSHEV_TRANSFORM_HPP
#include <cmath>
#include <typeinfo>
#include <type_traits>
#include <fftw3.h>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/chebyshev.hpp>

namespace boost { namespace math {

template<class Real>
class chebyshev_transform
{
public:
    template<class F>
    chebyshev_transform(const F& f, Real a, Real b, Real tol=500*std::numeric_limits<Real>::epsilon()): m_a(a), m_b(b)
    {
        static_assert(std::is_same<Real, float>::value || std::is_same<Real, double>::value
                   || std::is_same<Real, long double>::value,
                      "The Chebyshev transform only support float, double, and long double.\n");
        if (a >= b)
        {
            throw std::domain_error("a < b is required.\n");
        }
        using boost::math::constants::half;
        using boost::math::constants::pi;
        using std::cos;
        using std::abs;
        Real bma = (b-a)*half<Real>();
        Real bpa = (b+a)*half<Real>();
        size_t n = 4096;
        std::vector<Real> vf(n);
        m_coeffs.resize(n);


        // A C++17ism to make it so we can use different fftw_plans at compile time.
        // Is there a better way?
        if constexpr (std::is_same<Real, double>::value) {
            fftw_plan plan = fftw_plan_r2r_1d(n, vf.data(), m_coeffs.data(), FFTW_REDFT10, FFTW_ESTIMATE);
            Real inv_n = 1/static_cast<Real>(n);
            for(size_t j = 0; j < n; ++j)
            {
                Real y = cos(pi<Real>()*(j+half<Real>())/static_cast<Real>(n));
                vf[j] = f(y*bma + bpa)*inv_n;
            }

            fftw_execute_r2r(plan, vf.data(), m_coeffs.data());

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
            m_coeffs.resize(j+1);
        }
        // This is a way to get it to work generically on different precision types:
        /*else if constexpr (std::is_same<Real, float>::value) {
            fftwf_plan plan = fftwf_plan_r2r_1d(n, vf.data(), m_coeffs.data(), FFTW_REDFT10, FFTW_ESTIMATE);
            Real inv_n = 1/static_cast<Real>(n);
            for(size_t j = 0; j < n; ++j)
            {
                Real y = cos(pi<Real>()*(j+half<Real>())/static_cast<Real>(n));
                vf[j] = f(y*bma + bpa)*inv_n;
            }

            fftwf_execute_r2r(plan, vf.data(), m_coeffs.data());

            Real max_coeff = 0;
            for (auto const & coeff : m_coeffs)
            {
                if (abs(coeff) > max_coeff)
                {
                    max_coeff = abs(coeff);
                }
            }
            size_t j = m_coeffs.size() - 1;
            while (abs(m_coeffs[j])/max_coeff < 10*std::numeric_limits<Real>::epsilon())
            {
                --j;
            }
        }*/

    }

    Real operator()(Real x) const
    {
        using boost::math::constants::half;
        if (x > m_b || x < m_a)
        {
            throw std::domain_error("x not in [a, b]\n");
        }

        Real z = (2*x - m_a - m_b)/(m_b - m_a);
        return chebyshev_clenshaw_recurrence(m_coeffs.data(), m_coeffs.size(), z);
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

    void print_coefficients() const
    {
        std::cout << "{";
        for(auto const & coeff : m_coeffs) {
          std::cout << coeff << ", ";
        }
        std::cout << "}\n";
    }


private:
    std::vector<Real> m_coeffs;
    Real m_a;
    Real m_b;
};

}}
#endif
