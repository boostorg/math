//  (C) Copyright Nick Thompson 2018.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_DIFFERENTIATION_LANCZOS_SMOOTHING_HPP
#define BOOST_MATH_DIFFERENTIATION_LANCZOS_SMOOTHING_HPP
#include <vector>
#include <boost/assert.hpp>

namespace boost {
namespace math {
namespace differentiation {

namespace detail {
template <typename Real>
class discrete_legendre {
  public:
    discrete_legendre(int n) : m_n{n}
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
        m_qrm2 = 1;
        m_qrm1 = x;
        m_qrm2p = 0;
        m_qrm1p = 1;

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
        Real s =
            (m_r - 1) * (N * N - (m_r - 1) * (m_r - 1)) / Real(4 * m_n * m_n);
        Real tmp1 = ((2 * m_r - 1) * m_x * m_qrm1 - s * m_qrm2) / m_r;
        Real tmp2 = ((2 * m_r - 1) * (m_qrm1 + m_x * m_qrm1p) - s * m_qrm2p) / m_r;
        m_qrm2 = m_qrm1;
        m_qrm1 = tmp1;
        m_qrm2p = m_qrm1p;
        m_qrm1p = tmp2;
        ++m_r;
        return m_qrm1p;
    }


    Real operator()(Real x, int k)
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
        for (int r = 2; r <= k; ++r) {
            Real num = (r - 1) * (N * N - (r - 1) * (r - 1)) * qrm2;
            Real tmp = (2 * r - 1) * x * qrm1 - num / Real(4 * m_n * m_n);
            qrm2 = qrm1;
            qrm1 = tmp / r;
        }
        return qrm1;
    }

    Real prime(Real x, int k) {
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
        for (int r = 2; r <= k; ++r) {
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
    int m_n;
    Real m_qrm2;
    Real m_qrm1;
    Real m_qrm2p;
    Real m_qrm1p;
    int m_r;
    Real m_x;
};

template <class Real>
std::vector<Real> interior_filter(int n, int p) {
    // We could make the filter length n, as f[0] = 0,
    // but that'd make the indexing awkward when applying the filter.
    std::vector<Real> f(n + 1, 0);
    auto dlp = discrete_legendre<Real>(n);
    std::vector<Real> coeffs(p+1, std::numeric_limits<Real>::quiet_NaN());
    dlp.initialize_recursion(0);
    coeffs[1] = 1/dlp.norm_sq(1);
    for (int l = 3; l < p + 1; l += 2)
    {
        dlp.next_prime();
        coeffs[l] = dlp.next_prime()/ dlp.norm_sq(l);
    }

    for (size_t j = 1; j < f.size(); ++j)
    {
        Real arg = Real(j) / Real(n);
        dlp.initialize_recursion(arg);
        f[j] = coeffs[1]*arg;
        for (int l = 3; l <= p; l += 2)
        {
            dlp.next();
            f[j] += coeffs[l]*dlp.next();
        }
        f[j] /= (n * n);
    }
    return f;
}

template <class Real>
std::vector<Real> boundary_filter(int n, int p, int s) {
    std::vector<Real> f(2 * n + 1, 0);
    auto dlp = discrete_legendre<Real>(n);
    Real sn = Real(s) / Real(n);
    std::vector<Real> coeffs(p+1, std::numeric_limits<Real>::quiet_NaN());
    dlp.initialize_recursion(sn);
    coeffs[0] = 0;
    coeffs[1] = 1/dlp.norm_sq(1);
    for (int l = 2; l < p + 1; ++l)
    {
        // Calculation of the norms is common to all filters,
        // so it seems like an obvious optimization target.
        // I tried this: The spent in computing the norms time is not negligible,
        // but still a small fraction of the total compute time.
        // Hence I'm not refactoring out these norm calculations.
        coeffs[l] = dlp.next_prime()/ dlp.norm_sq(l);
    }

    for (int k = 0; k < f.size(); ++k)
    {
        int j = k - n;
        f[k] = 0;
        Real arg = Real(j) / Real(n);
        dlp.initialize_recursion(arg);
        f[k] = coeffs[1]*arg;
        for (int l = 2; l <= p; ++l)
        {
            f[k] += coeffs[l]*dlp.next();
        }
        f[k] /= (n * n);
    }
    return f;
}

} // namespace detail

template <class RandomAccessContainer>
class lanczos_derivative {
public:
    using Real = typename RandomAccessContainer::value_type;
    lanczos_derivative(RandomAccessContainer const &v,
                       Real spacing = 1,
                       int filter_length = 18,
                       int approximation_order = 3)
        : m_v{v}, dt{spacing}
    {
        BOOST_ASSERT_MSG(approximation_order <= 2 * filter_length,
                         "The approximation order must be <= 2n");
        BOOST_ASSERT_MSG(approximation_order >= 2,
                         "The approximation order must be >= 2");
        BOOST_ASSERT_MSG(spacing > 0, "Spacing between samples must be > 0.");
        using std::size;
        BOOST_ASSERT_MSG(size(v) >= filter_length,
            "Vector must be at least as long as the filter length");
        m_f = detail::interior_filter<Real>(filter_length, approximation_order);

        boundary_filters.resize(filter_length);
        for (size_t i = 0; i < filter_length; ++i)
        {
            // s = i - n;
            boundary_filters[i] = detail::boundary_filter<Real>(
                filter_length, approximation_order, i - filter_length);
        }
    }

    void reset_data(RandomAccessContainer const &v)
    {
        using std::size;
        BOOST_ASSERT_MSG(size(v) >= m_f.size(), "Vector must be at least as long as the filter length");
        m_v = v;
    }

    void reset_spacing(Real spacing)
    {
        BOOST_ASSERT_MSG(spacing > 0, "Spacing between samples must be > 0.");
        dt = spacing;
    }

    Real spacing() const
    {
        return dt;
    }

    Real operator[](size_t i) const
    {
        using std::size;
        if (i >= m_f.size() - 1 && i <= size(m_v) - m_f.size())
        {
            Real dv = 0;
            for (size_t j = 1; j < m_f.size(); ++j)
            {
                dv += m_f[j] * (m_v[i + j] - m_v[i - j]);
            }
            return dv / dt;
        }

        // m_f.size() = N+1
        if (i < m_f.size() - 1)
        {
            auto &f = boundary_filters[i];
            Real dv = 0;
            for (size_t j = 0; j < f.size(); ++j) {
                dv += f[j] * m_v[j];
            }
            return dv / dt;
        }

        if (i > size(m_v) - m_f.size() && i < size(m_v))
        {
            int k = size(m_v) - 1 - i;
            auto &f = boundary_filters[k];
            Real dv = 0;
            for (size_t j = 0; j < f.size(); ++j)
            {
                dv += f[j] * m_v[m_v.size() - 1 - j];
            }
            return -dv / dt;
        }

        // OOB access:
        BOOST_ASSERT_MSG(false, "Out of bounds access in Lanczos derivative");
        return std::numeric_limits<Real>::quiet_NaN();
    }

private:
    const RandomAccessContainer &m_v;
    std::vector<Real> m_f;
    std::vector<std::vector<Real>> boundary_filters;
    Real dt;
};

// We can also implement lanczos_acceleration, but let's get the API for lanczos_derivative nailed down before doing so.

} // namespace differentiation
} // namespace math
} // namespace boost
#endif
