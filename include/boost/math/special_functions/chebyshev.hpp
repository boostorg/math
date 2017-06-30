//  (C) Copyright Nick Thompson 2017.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_SPECIAL_CHEBYSHEV_HPP
#define BOOST_MATH_SPECIAL_CHEBYSHEV_HPP
#include <cmath>
#include <boost/math/constants/constants.hpp>

namespace boost { namespace math {

template<class Real>
inline Real chebyshev_next(Real const & x, Real const & Tn, Real const & Tn_1)
{
    return 2*x*Tn - Tn_1;
}

namespace detail {

template<class Real>
std::vector<Real> discrete_cosine_transform_type_II(const Real* const x, size_t n)
{
    using std::cos;
    using boost::math::constants::pi;
    using boost::math::constants::half;
    std::vector<Real> X(n, 0);
    Real scale = static_cast<Real>(2)/static_cast<Real>(n);
    for(size_t k = 0; k < n; ++k)
    {
        Real g = pi<Real>()*k/n;
        for (size_t j = 0; j < n; ++j)
        {
           X[k] += x[k]*cos(g*(j+half<Real>()));
        }
        X[k] *= scale;
    }
    return X;
}

template<class Real, bool second=false>
inline Real chebyshev_imp(unsigned n, Real const & x)
{
    using std::cosh;
    using std::acosh;
    using std::pow;
    Real T0 = 1;
    Real T1;
    if (second)
    {
        if (x > 1 || x < -1)
        {
            Real t = sqrt(x*x -1);
            return (pow(x+t, n+1) - pow(x-t, n+1))/(2*t);
        }
        T1 = 2*x;
    }
    else
    {
        if (x > 1)
        {
            return cosh(n*acosh(x));
        }
        if (x < -1)
        {
            if (n & 1)
            {
                return -cosh(n*acosh(-x));
            }
            else
            {
                return cosh(n*acosh(-x));
            }
        }
        T1 = x;
    }

    if (n == 0)
    {
        return T0;
    }

    unsigned l = 1;
    while(l < n)
    {
       std::swap(T0, T1);
       T1 = boost::math::chebyshev_next(x, T0, T1);
       ++l;
    }
    return T1;
}
} // namespace detail

template<class Real>
Real chebyshev_t(unsigned n, Real const & x)
{
    return detail::chebyshev_imp<Real, false>(n, x);
}

template<class Real>
Real chebyshev_u(unsigned n, Real const & x)
{
    return detail::chebyshev_imp<Real, true>(n, x);
}

template<class Real>
Real chebyshev_t_prime(unsigned n, Real const & x)
{
    if (n == 0)
    {
        return 0;
    }
    return n*detail::chebyshev_imp<Real, true>(n - 1, x);
}

/*
 * This is Algorithm 3.1 of
 * Gil, Amparo, Javier Segura, and Nico M. Temme.
 * Numerical methods for special functions.
 * Society for Industrial and Applied Mathematics, 2007.
 * https://www.siam.org/books/ot99/OT99SampleChapter.pdf
 * However, our definition of c0 differs . . .
 */
template<class Real>
inline Real chebyshev_clenshaw_recurrence(const Real* const c, size_t length, Real x)
{
    using boost::math::constants::half;
    Real b2 = 0;
    Real b1 = c[length -1];
    for(size_t j = length - 2; j >= 1; --j)
    {
        Real tmp = 2*x*b1 - b2 + c[j];
        b2 = b1;
        b1 = tmp;
    }
    return x*b1 - b2 + half<Real>()*c[0];
}

template<class Real>
class chebyshev_transform
{
public:
    template<class F>
    chebyshev_transform(const F& f, Real a, Real b): m_a(a), m_b(b)
    {
        if (a >= b)
        {
            throw std::domain_error("a < b is required.\n");
        }
        using boost::math::constants::half;
        using boost::math::constants::pi;
        using std::cos;
        using std::abs;
        size_t n = 1024;
        std::vector<Real> vf(n);
        Real bma = (b-a)*half<Real>();
        Real bpa = (b+a)*half<Real>();
        for(size_t k = 0; k < n; ++k)
        {
            Real y = cos(pi<Real>()*(k+half<Real>())/static_cast<Real>(n));
            vf[k] = f(y*bma + bpa);
        }
        m_coeffs = detail::discrete_cosine_transform_type_II(vf.data(), vf.size());
        /*Real scale = static_cast<Real>(2)/static_cast<Real>(n);
        m_coeffs.resize(n);
        Real max_coeff = 0;
        for(size_t j = 0; j < n; ++j)
        {
            Real compensator = 0;
            Real sum = 0;
            for(size_t k = 0; k < n; ++k)
            {
                Real y = vf[k]*cos(pi<Real>()*j*(k+half<Real>())/n);
                Real t = sum + y;
                compensator = (t - sum) - y;
                sum = t;
            }
            m_coeffs[j] = scale*sum;
            if (abs(m_coeffs[j]) > max_coeff)
            {
                max_coeff = abs(m_coeffs[j]);
            }
        }
        for(auto c : m_coeffs)
        {
           std::cout << c << "\t";
        }
        std::cout << std::endl;
        size_t j = m_coeffs.size() - 1;
        while (abs(m_coeffs[j])/max_coeff < 10*std::numeric_limits<Real>::epsilon())
        {
            --j;
        }

        n = j + 2;
        vf.resize(n);
        for(size_t k = 0; k < n; ++k)
        {
            Real y = cos(pi<Real>()*(k+half<Real>())/static_cast<Real>(n));
            vf[k] = f(y*bma + bpa);
        }
        scale = static_cast<Real>(2)/static_cast<Real>(n);
        m_coeffs.resize(n);
        for(size_t j = 0; j < n; ++j)
        {
            Real compensator = 0;
            Real sum = 0;
            for(size_t k = 0; k < n; ++k)
            {
                Real y = vf[k]*cos(pi<Real>()*j*(k+half<Real>())/n);
                Real t = sum + y;
                compensator = (t - sum) - y;
                sum = t;
            }
            m_coeffs[j] = scale*sum;
        }
        std::cout << "Needed " << m_coeffs.size() << " coefficients\n";
        for(auto c : m_coeffs)
        {
           std::cout << c << "\t";
        }
        std::cout << std::endl;*/

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

    // Integral from a <= x1 <= x2 <= b
    Real integrate(Real x1, Real x2) const
    {
        if (x1 < m_a || x2 > m_b)
        {
            throw std::domain_error("Cannot integrate outside boundary of interpolation.\n");
        }
        return std::numeric_limits<Real>::quiet_Nan();
    }

    const std::vector<Real>& coefficients() const
    {
        return m_coeffs;
    }

    Real prime(Real x) const
    {
        // T_{n+1} = 2*x*T_{n} - T_{n-1}
        using boost::math::constants::half;
        Real yp = 0;
        Real z = (2*x - m_a - m_b)/(m_b - m_a);
        for (size_t i = 1; i < m_coeffs.size(); ++i)
        {
            yp += m_coeffs[i]*chebyshev_t_prime(i, z);
        }
        return 2*yp/(m_b - m_a);
    }


private:
    std::vector<Real> m_coeffs;
    Real m_a;
    Real m_b;
};

}}
#endif
