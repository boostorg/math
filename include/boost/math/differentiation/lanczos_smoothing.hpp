//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DIFFERENTIATION_LANCZOS_SMOOTHING_HPP
#define BOOST_MATH_DIFFERENTIATION_LANCZOS_SMOOTHING_HPP
#include <vector>
#include <boost/assert.hpp>

namespace boost::math::differentiation {

namespace detail {
template <typename Real>
class discrete_legendre {
  public:
    explicit discrete_legendre(size_t n) : m_n{n}, m_r{2},
                                           m_x{std::numeric_limits<Real>::quiet_NaN()},
                                           m_qrm2{std::numeric_limits<Real>::quiet_NaN()},
                                           m_qrm1{std::numeric_limits<Real>::quiet_NaN()},
                                           m_qrm2p{std::numeric_limits<Real>::quiet_NaN()},
                                           m_qrm1p{std::numeric_limits<Real>::quiet_NaN()},
                                           m_qrm2pp{std::numeric_limits<Real>::quiet_NaN()},
                                           m_qrm1pp{std::numeric_limits<Real>::quiet_NaN()}
    {
        // The integer n indexes a family of discrete Legendre polynomials indexed by k <= 2*n
    }

    Real norm_sq(int r)
    {
        Real prod = Real(2) / Real(2 * r + 1);
        for (int k = -r; k <= r; ++k) {
            prod *= Real(2 * m_n + 1 + k) / Real(2 * m_n);
        }
        return prod;
    }

    void initialize_recursion(Real x)
    {
        using std::abs;
        BOOST_ASSERT_MSG(abs(x) <= 1, "Three term recurrence is stable only for |x| <=1.");
        m_qrm2 = 1;
        m_qrm1 = x;
        // Derivatives:
        m_qrm2p = 0;
        m_qrm1p = 1;
        // Second derivatives:
        m_qrm2pp = 0;
        m_qrm1pp = 0;

        m_r = 2;
        m_x = x;
    }

    Real next()
    {
        Real N = 2 * m_n + 1;
        Real num = (m_r - 1) * (N * N - (m_r - 1) * (m_r - 1)) * m_qrm2;
        Real tmp = (2 * m_r - 1) * m_x * m_qrm1 - num / Real(4 * m_n * m_n);
        m_qrm2 = m_qrm1;
        m_qrm1 = tmp / m_r;
        ++m_r;
        return m_qrm1;
    }

    Real next_prime()
    {
        Real N = 2 * m_n + 1;
        Real s = (m_r - 1) * (N * N - (m_r - 1) * (m_r - 1)) / Real(4 * m_n * m_n);
        Real tmp1 = ((2 * m_r - 1) * m_x * m_qrm1 - s * m_qrm2) / m_r;
        Real tmp2 = ((2 * m_r - 1) * (m_qrm1 + m_x * m_qrm1p) - s * m_qrm2p) / m_r;
        m_qrm2 = m_qrm1;
        m_qrm1 = tmp1;
        m_qrm2p = m_qrm1p;
        m_qrm1p = tmp2;
        ++m_r;
        return m_qrm1p;
    }

    Real next_dbl_prime()
    {
        Real N = 2*m_n + 1;
        Real trm1 = 2*m_r - 1;
        Real s = (m_r - 1) * (N * N - (m_r - 1) * (m_r - 1)) / Real(4 * m_n * m_n);
        Real rqrpp = 2*trm1*m_qrm1p + trm1*m_x*m_qrm1pp - s*m_qrm2pp;
        Real tmp1 = ((2 * m_r - 1) * m_x * m_qrm1 - s * m_qrm2) / m_r;
        Real tmp2 = ((2 * m_r - 1) * (m_qrm1 + m_x * m_qrm1p) - s * m_qrm2p) / m_r;
        m_qrm2 = m_qrm1;
        m_qrm1 = tmp1;
        m_qrm2p = m_qrm1p;
        m_qrm1p = tmp2;
        m_qrm2pp = m_qrm1pp;
        m_qrm1pp = rqrpp/m_r;
        ++m_r;
        return m_qrm1pp;
    }

    Real operator()(Real x, size_t k)
    {
        BOOST_ASSERT_MSG(k <= 2 * m_n, "r <= 2n is required.");
        if (k == 0)
        {
            return 1;
        }
        if (k == 1)
        {
            return x;
        }
        Real qrm2 = 1;
        Real qrm1 = x;
        Real N = 2 * m_n + 1;
        for (size_t r = 2; r <= k; ++r) {
            Real num = (r - 1) * (N * N - (r - 1) * (r - 1)) * qrm2;
            Real tmp = (2 * r - 1) * x * qrm1 - num / Real(4 * m_n * m_n);
            qrm2 = qrm1;
            qrm1 = tmp / r;
        }
        return qrm1;
    }

    Real prime(Real x, size_t k) {
        BOOST_ASSERT_MSG(k <= 2 * m_n, "r <= 2n is required.");
        if (k == 0) {
            return 0;
        }
        if (k == 1) {
            return 1;
        }
        Real qrm2 = 1;
        Real qrm1 = x;
        Real qrm2p = 0;
        Real qrm1p = 1;
        Real N = 2 * m_n + 1;
        for (size_t r = 2; r <= k; ++r) {
            Real s =
                (r - 1) * (N * N - (r - 1) * (r - 1)) / Real(4 * m_n * m_n);
            Real tmp1 = ((2 * r - 1) * x * qrm1 - s * qrm2) / r;
            Real tmp2 = ((2 * r - 1) * (qrm1 + x * qrm1p) - s * qrm2p) / r;
            qrm2 = qrm1;
            qrm1 = tmp1;
            qrm2p = qrm1p;
            qrm1p = tmp2;
        }
        return qrm1p;
    }

  private:
    size_t m_n;
    size_t m_r;
    Real m_x;
    Real m_qrm2;
    Real m_qrm1;
    Real m_qrm2p;
    Real m_qrm1p;
    Real m_qrm2pp;
    Real m_qrm1pp;
};

template <class Real>
std::vector<Real> interior_filter(size_t n, size_t p) {
    // We could make the filter length n, as f[0] = 0,
    // but that'd make the indexing awkward when applying the filter.
    std::vector<Real> f(n + 1, 0);
    auto dlp = discrete_legendre<Real>(n);
    std::vector<Real> coeffs(p+1, std::numeric_limits<Real>::quiet_NaN());
    dlp.initialize_recursion(0);
    coeffs[1] = 1/dlp.norm_sq(1);
    for (size_t l = 3; l < p + 1; l += 2)
    {
        dlp.next_prime();
        coeffs[l] = dlp.next_prime()/ dlp.norm_sq(l);
    }

    for (size_t j = 1; j < f.size(); ++j)
    {
        Real arg = Real(j) / Real(n);
        dlp.initialize_recursion(arg);
        f[j] = coeffs[1]*arg;
        for (size_t l = 3; l <= p; l += 2)
        {
            dlp.next();
            f[j] += coeffs[l]*dlp.next();
        }
        f[j] /= (n * n);
    }
    return f;
}

template <class Real>
std::vector<Real> boundary_filter(size_t n, size_t p, int64_t s)
{
    std::vector<Real> f(2 * n + 1, 0);
    auto dlp = discrete_legendre<Real>(n);
    Real sn = Real(s) / Real(n);
    std::vector<Real> coeffs(p+1, std::numeric_limits<Real>::quiet_NaN());
    dlp.initialize_recursion(sn);
    coeffs[0] = 0;
    coeffs[1] = 1/dlp.norm_sq(1);
    for (size_t l = 2; l < p + 1; ++l)
    {
        // Calculation of the norms is common to all filters,
        // so it seems like an obvious optimization target.
        // I tried this: The spent in computing the norms time is not negligible,
        // but still a small fraction of the total compute time.
        // Hence I'm not refactoring out these norm calculations.
        coeffs[l] = dlp.next_prime()/ dlp.norm_sq(l);
    }

    for (size_t k = 0; k < f.size(); ++k)
    {
        Real j = Real(k) - Real(n);
        f[k] = 0;
        Real arg = j/Real(n);
        dlp.initialize_recursion(arg);
        f[k] = coeffs[1]*arg;
        for (size_t l = 2; l <= p; ++l)
        {
            f[k] += coeffs[l]*dlp.next();
        }
        f[k] /= (n * n);
    }
    return f;
}

template <class Real>
std::vector<Real> acceleration_boundary_filter(size_t n, size_t p, int64_t s)
{
    BOOST_ASSERT_MSG(p <= 2*n, "Approximation order must be <= 2*n");
    BOOST_ASSERT_MSG(p > 2, "Approximation order must be > 2");
    std::vector<Real> f(2 * n + 1, 0);
    auto dlp = discrete_legendre<Real>(n);
    Real sn = Real(s) / Real(n);
    std::vector<Real> coeffs(p+2, std::numeric_limits<Real>::quiet_NaN());
    dlp.initialize_recursion(sn);
    coeffs[0] = 0;
    coeffs[1] = 0;
    for (size_t l = 2; l < p + 2; ++l)
    {
        coeffs[l] = dlp.next_dbl_prime()/ dlp.norm_sq(l);
    }

    for (size_t k = 0; k < f.size(); ++k)
    {
        Real j = Real(k) - Real(n);
        f[k] = 0;
        Real arg = j/Real(n);
        dlp.initialize_recursion(arg);
        f[k] = coeffs[1]*arg;
        for (size_t l = 2; l <= p; ++l)
        {
            f[k] += coeffs[l]*dlp.next();
        }
        f[k] /= (n * n * n);
    }
    return f;
}


} // namespace detail

template <typename Real, size_t order = 1>
class discrete_lanczos_derivative {
public:
    discrete_lanczos_derivative(Real const & spacing,
                                size_t n = 18,
                                size_t approximation_order = 3)
        : m_dt{spacing}
    {
        static_assert(!std::is_integral_v<Real>, "Spacing must be a floating point type.");
        BOOST_ASSERT_MSG(spacing > 0, "Spacing between samples must be > 0.");

        if constexpr (order == 1)
        {
            BOOST_ASSERT_MSG(approximation_order <= 2 * n,
                             "The approximation order must be <= 2n");
            BOOST_ASSERT_MSG(approximation_order >= 2,
                             "The approximation order must be >= 2");
            m_f = detail::interior_filter<Real>(n, approximation_order);

            m_boundary_filters.resize(n);
            for (size_t i = 0; i < n; ++i)
            {
                // s = i - n;
                int64_t s = static_cast<int64_t>(i) - static_cast<int64_t>(n);
                m_boundary_filters[i] = detail::boundary_filter<Real>(n, approximation_order, s);
            }
        }
        else if constexpr (order == 2)
        {
            auto f = detail::acceleration_boundary_filter<Real>(n, approximation_order, 0);
            m_f.resize(n+1);
            for (size_t i = 0; i < m_f.size(); ++i)
            {
                m_f[i] = f[i+n];
            }
            m_boundary_filters.resize(n);
            for (size_t i = 0; i < n; ++i)
            {
                int64_t s = static_cast<int64_t>(i) - static_cast<int64_t>(n);
                m_boundary_filters[i] = detail::acceleration_boundary_filter<Real>(n, approximation_order, s);
            }
        }
        else
        {
            BOOST_ASSERT_MSG(false, "Derivatives of order 3 and higher are not implemented.");
        }
    }

    void reset_spacing(Real const & spacing)
    {
        BOOST_ASSERT_MSG(spacing > 0, "Spacing between samples must be > 0.");
        m_dt = spacing;
    }

    Real spacing() const
    {
        return m_dt;
    }

    template<class RandomAccessContainer>
    Real operator()(RandomAccessContainer const & v, size_t i) const
    {
        static_assert(std::is_same_v<typename RandomAccessContainer::value_type, Real>,
                      "The type of the values in the vector provided does not match the type in the filters.");
        using std::size;
        BOOST_ASSERT_MSG(size(v) >= m_boundary_filters[0].size(),
            "Vector must be at least as long as the filter length");

        if constexpr (order==1)
        {
            if (i >= m_f.size() - 1 && i <= size(v) - m_f.size())
            {
                Real dv = 0;
                for (size_t j = 1; j < m_f.size(); ++j)
                {
                    dv += m_f[j] * (v[i + j] - v[i - j]);
                }
                return dv / m_dt;
            }

            // m_f.size() = N+1
            if (i < m_f.size() - 1)
            {
                auto &bf = m_boundary_filters[i];
                Real dv = 0;
                for (size_t j = 0; j < bf.size(); ++j)
                {
                    dv += bf[j] * v[j];
                }
                return dv / m_dt;
            }

            if (i > size(v) - m_f.size() && i < size(v))
            {
                int k = size(v) - 1 - i;
                auto &bf = m_boundary_filters[k];
                Real dv = 0;
                for (size_t j = 0; j < bf.size(); ++j)
                {
                    dv += bf[j] * v[size(v) - 1 - j];
                }
                return -dv / m_dt;
            }
        }
        else if constexpr (order==2)
        {
            if (i >= m_f.size() - 1 && i <= size(v) - m_f.size())
            {
                Real d2v = m_f[0]*v[i];
                for (size_t j = 1; j < m_f.size(); ++j)
                {
                    d2v += m_f[j] * (v[i + j] + v[i - j]);
                }
                return d2v / (m_dt*m_dt);
            }

            // m_f.size() = N+1
            if (i < m_f.size() - 1)
            {
                auto &bf = m_boundary_filters[i];
                Real d2v = 0;
                for (size_t j = 0; j < bf.size(); ++j)
                {
                    d2v += bf[j] * v[j];
                }
                return d2v / (m_dt*m_dt);
            }

            if (i > size(v) - m_f.size() && i < size(v))
            {
                int k = size(v) - 1 - i;
                auto &bf = m_boundary_filters[k];
                Real d2v = 0;
                for (size_t j = 0; j < bf.size(); ++j)
                {
                    d2v += bf[j] * v[size(v) - 1 - j];
                }
                return d2v / (m_dt*m_dt);
            }
        }

        // OOB access:
        BOOST_ASSERT_MSG(false, "Out of bounds access in Lanczos derivative");
        return std::numeric_limits<Real>::quiet_NaN();
    }

    template<class RandomAccessContainer>
    RandomAccessContainer operator()(RandomAccessContainer const & v) const
    {
        static_assert(std::is_same_v<typename RandomAccessContainer::value_type, Real>,
                      "The type of the values in the vector provided does not match the type in the filters.");
        using std::size;
        BOOST_ASSERT_MSG(size(v) >= m_boundary_filters[0].size(),
            "Vector must be at least as long as the filter length");

        RandomAccessContainer w(size(v));
        if constexpr (order==1)
        {
            for (size_t i = 0; i < m_f.size() - 1; ++i)
            {
                auto &bf = m_boundary_filters[i];
                Real dv = 0;
                for (size_t j = 0; j < bf.size(); ++j)
                {
                    dv += bf[j] * v[j];
                }
                w[i] = dv / m_dt;
            }

            for(size_t i = m_f.size() - 1; i <= size(v) - m_f.size(); ++i)
            {
                Real dv = 0;
                for (size_t j = 1; j < m_f.size(); ++j)
                {
                    dv += m_f[j] * (v[i + j] - v[i - j]);
                }
                w[i] = dv / m_dt;
            }


            for(size_t i = size(v) - m_f.size() + 1; i < size(v); ++i)
            {
                int k = size(v) - 1 - i;
                auto &f = m_boundary_filters[k];
                Real dv = 0;
                for (size_t j = 0; j < f.size(); ++j)
                {
                    dv += f[j] * v[size(v) - 1 - j];
                }
                w[i] = -dv / m_dt;
            }
        }
        else if constexpr (order==2)
        {
            // m_f.size() = N+1
            for (size_t i = 0; i < m_f.size() - 1; ++i)
            {
                auto &bf = m_boundary_filters[i];
                Real d2v = 0;
                for (size_t j = 0; j < bf.size(); ++j)
                {
                    d2v += bf[j] * v[j];
                }
                w[i] = d2v / (m_dt*m_dt);
            }

            for (size_t i = m_f.size() - 1; i <= size(v) - m_f.size(); ++i)
            {
                Real d2v = m_f[0]*v[i];
                for (size_t j = 1; j < m_f.size(); ++j)
                {
                    d2v += m_f[j] * (v[i + j] + v[i - j]);
                }
                w[i] = d2v / (m_dt*m_dt);
            }


            for (size_t i = size(v) - m_f.size() + 1; i < size(v); ++i)
            {
                int k = size(v) - 1 - i;
                auto &bf = m_boundary_filters[k];
                Real d2v = 0;
                for (size_t j = 0; j < bf.size(); ++j)
                {
                    d2v += bf[j] * v[size(v) - 1 - j];
                }
                w[i] = d2v / (m_dt*m_dt);
            }
        }

        return w;
    }

private:
    std::vector<Real> m_f;
    std::vector<std::vector<Real>> m_boundary_filters;
    Real m_dt;
};

} // namespaces
#endif
