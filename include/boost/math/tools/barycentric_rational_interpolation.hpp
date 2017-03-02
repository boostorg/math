/*
 *  Copyright Nick Thompson, 2017
 *  Use, modification and distribution are subject to the
 *  Boost Software License, Version 1.0. (See accompanying file
 *  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
 *
 *  Given N samples (t_i, y_i) which are irregularly spaced, this routine constructs an
 *  interpolant s which is constructed in O(N) time, occupies O(N) space, and can be evaluated in O(N) time.
 *  The interpolation is stable, unless one point is incredibly close to another, and the next point is incredibly far.
 *  The measure of this stability is the "local mesh ratio", which can be queried from the routine.
 *  Pictorially, the following t_i spacing is bad (has a high local mesh ratio)
 *  ||             |      | |                           |
 *  and this t_i spacing is good (has a low local mesh ratio)
 *  |   |      |    |     |    |        |    |  |    |
 *
 *
 *  If f is C^{d+2}, then the interpolant is O(h^(d+1)) accurate, where d is the interpolation order.
 *  A disadvantage of this interpolant is that it does not reproduce rational functions; for example, 1/(1+x^2) is not interpolated exactly.
 *
 *  This routine follows "Barycentric Rational Approximation with no Poles and High Rate of Approximation"
 *  by Michael S. Floater, Kai Hormann
 *  http://www.mn.uio.no/math/english/people/aca/michaelf/papers/rational.pdf
 *  Official Citation:
 *  Floater, M.S. & Hormann, K. Numer. Math. (2007) 107: 315. doi:10.1007/s00211-007-0093-y
 *
 *  Floater's algorithm has been discussed in detail in Numerical Recipes, 3rd edition, section 3.4.1.
 */

#ifndef BOOST_MATH_TOOLS_BARYCENTRIC_RATIONAL_INTERPOLATION_HPP
#define BOOST_MATH_TOOLS_BARYCENTRIC_RATIONAL_INTERPOLATION_HPP

#include <vector>
#include <boost/type_index.hpp>
#include <boost/format.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

namespace boost{ namespace math{ namespace tools{

template<class Real>
class barycentric_rational
{
public:
    barycentric_rational(const Real* const x, const Real* const y, size_t n, size_t approximation_order = 3);

    Real operator()(Real t) const;

    // The barycentric weights are not really that interesting; except to the unit tests!
    Real weight(size_t i) const { return m_w[i]; }

private:
    // Technically, we do not need to copy over y to m_y, or x to m_x.
    // We could simply store a pointer. However, in doing so,
    // we'd need to make sure the user managed the lifetime of m_y,
    // and didn't alter its data. Since we are unlikely to run out of
    // memory for a linearly scaling algorithm, it seems best just to make a copy.
    std::vector<Real> m_y;
    std::vector<Real> m_x;
    std::vector<Real> m_w;
};

template<class Real>
barycentric_rational<Real>::barycentric_rational(const Real* const x, const Real* const y, size_t n, size_t approximation_order)
{
    using std::abs;

    if (approximation_order >= n)
    {
        throw std::domain_error("Approximation order must be < data length.");
    }

    if (x == nullptr)
    {
        throw std::domain_error("Independent variable passed to barycentric_rational is a nullptr");
    }

    if (y == nullptr)
    {
        throw std::domain_error("Dependent variable passed to barycentric_rational is a nullptr");
    }

    // Big sad memcpy to make sure the object is easy to use.
    m_x.resize(n);
    m_y.resize(n);
    for(size_t i = 0; i < n; ++i)
    {
        // But if we're going to do a memcpy, we can do some error checking which is inexpensive relative to the copy:
        if(boost::math::isnan(x[i]))
        {
            boost::format fmtr = boost::format("x[%1%] is a NAN") % i;
            throw std::domain_error(fmtr.str());
        }

        if(boost::math::isnan(y[i]))
        {
            boost::format fmtr = boost::format("y[%1%] is a NAN") % i;
            throw std::domain_error(fmtr.str());
        }

        m_x[i] = x[i];
        m_y[i] = y[i];
    }

    m_w.resize(n, 0);
    for(int64_t k = 0; k < n; ++k)
    {
        int64_t i_min = std::max(k - (int64_t) approximation_order, (int64_t) 0);
        int64_t i_max = k;
        if (k >= n - approximation_order)
        {
            i_max = n - approximation_order - 1;
        }

        for(int64_t i = i_min; i <= i_max; ++i)
        {
            Real inv_product = 1;
            int64_t j_max = std::min(i + approximation_order, n - 1);
            for(int64_t j = i; j <= j_max; ++j)
            {
                if (j == k)
                {
                    continue;
                }

                Real diff = m_x[k] - m_x[j];
                if (abs(diff) < std::numeric_limits<Real>::epsilon())
                {
                    boost::format fmtr = boost::format("Spacing between  x[%1%] and x[%2%] is %3%, which is smaller than the epsilon of %4%") % k % i % diff % boost::typeindex::type_id<Real>().pretty_name();
                    std::cout << typeid(fmtr).name() << std::endl;
                    throw std::logic_error(fmtr.str());
                }
                inv_product *= diff;
            }
            if (i % 2 == 0)
            {
                m_w[k] += 1/inv_product;
            }
            else
            {
                m_w[k] -= 1/inv_product;
            }
        }
    }
}

template<class Real>
Real barycentric_rational<Real>::operator()(Real x) const
{
    Real numerator = 0;
    Real denominator = 0;
    for(size_t i = 0; i < m_x.size(); ++i)
    {
        // Presumably we should see if the accuracy is improved by using ULP distance of say, 5 here, instead of testing for floating point equality.
        // However, it has been shown that if x approx x_i, but x != x_i, then inaccuracy in the numerator cancels the inaccuracy in the denominator,
        // and the result is fairly accurate. See: http://epubs.siam.org/doi/pdf/10.1137/S0036144502417715
        if (x == m_x[i])
        {
            return m_y[i];
        }
        Real t = m_w[i]/(x - m_x[i]);
        numerator += t*m_y[i];
        denominator += t;
    }
    return numerator/denominator;
}

/*
 * A formula for computing the derivative of the barycentric representation is given in
 * "Some New Aspects of Rational Interpolation", by Claus Schneider and Wilhelm Werner,
 * Mathematics of Computation, v47, number 175, 1986.
 * http://www.ams.org/journals/mcom/1986-47-175/S0025-5718-1986-0842136-8/S0025-5718-1986-0842136-8.pdf
 * However, this requires a lot of machinery which is not built into the library at present.
 * So we wait until there is a requirement to interpolate the derivative.
 */

}}}
#endif
