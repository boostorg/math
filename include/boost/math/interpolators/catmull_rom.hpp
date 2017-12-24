// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// This computes the Catmull-Rom spline from a list of points.

#ifndef BOOST_MATH_INTERPOLATORS_CATMULL_ROM
#define BOOST_MATH_INTERPOLATORS_CATMULL_ROM

#include <cmath>
#include <cassert>
#include <algorithm>

namespace boost{ namespace math{

    namespace detail
    {
        template<class Real, class Point, size_t dimension>
        Real alpha_distance(Point const & p1, Point const & p2, Real alpha)
        {
            using std::pow;
            Real dsq = 0;
            for (size_t i = 0; i < dimension; ++i)
            {
                Real dx = p1[i] - p2[i];
                dsq += dx*dx;
            }
            return pow(dsq, alpha/2);
        }
    }

template <class Real, class Point, size_t dimension>
class catmull_rom
{
public:

    catmull_rom(const Point* const points, size_t num_pnts, bool closed = false, Real alpha = (Real) 1/ (Real) 2);

    Real max_parameter() const
    {
        return m_max_s;
    }

    Real parameter_at_point(size_t i) const
    {
        return m_s[i+1];
    }

    Point operator()(Real s) const;


    Point prime(Real s) const;

private:
    std::vector<Point> m_pnts;
    std::vector<Real> m_s;
    Real m_max_s;
};

template<class Real, class Point, size_t dimension>
catmull_rom<Real, Point, dimension>::catmull_rom(const Point* const points, size_t num_pnts, bool closed, Real alpha)
{
    static_assert(dimension > 0, "The dimension of the Catmull-Rom spline must be > 0\n");
    if (num_pnts < 4)
    {
        throw std::domain_error("The Catmull-Rom curve requires at least 4 points.\n");
    }
    using std::abs;
    m_s.resize(num_pnts+3);
    m_pnts.resize(num_pnts+3);

    m_pnts[0] = points[num_pnts-1];
    for (size_t i = 0; i < num_pnts; ++i)
    {
        m_pnts[i+1] = points[i];
    }
    m_pnts[num_pnts+1] = points[0];
    m_pnts[num_pnts+2] = points[1];
    m_s[0] = -detail::alpha_distance<Real, Point, dimension>(m_pnts[0], m_pnts[1], alpha);
    if (abs(m_s[0]) < std::numeric_limits<Real>::epsilon())
    {
        throw std::domain_error("The first and last point should not be the same.\n");
    }
    m_s[1] = 0;
    for (size_t i = 2; i < m_s.size(); ++i)
    {
        Real d = detail::alpha_distance<Real, Point, dimension>(m_pnts[i], m_pnts[i-1], alpha);
        if (abs(d) < std::numeric_limits<Real>::epsilon())
        {
            throw std::domain_error("The control points of the Catmull-Rom curve are too close together; this will lead to ill-conditioning.\n");
        }
        m_s[i] = m_s[i-1] + d;
    }
    if(closed)
    {
        m_max_s = m_s[num_pnts+1];
    }
    else
    {
        m_max_s = m_s[num_pnts];
    }
}


template<class Real, class Point, size_t dimension>
Point catmull_rom<Real, Point, dimension>::operator()(Real s) const
{
    if (s < 0 || s > m_max_s)
    {
        throw std::domain_error("Parameter outside bounds.\n");
    }
    auto it = std::upper_bound(m_s.begin(), m_s.end(), s);
    //Now *it >= s. We want the index such that m_s[i] <= s < m_s[i+1]:
    size_t i = std::distance(m_s.begin(), it - 1);
    // We'll keep the assert in here a while until we feel good that we've understood this algorithm.
    //assert(m_s[i] <= s && s < m_s[i+1]);

    // Only denom21 is used twice:
    Real denom21 = 1/(m_s[i+1] - m_s[i]);
    Real s0s = m_s[i-1] - s;
    Real s1s = m_s[i] - s;
    Real s2s = m_s[i+1] - s;
    Real s3s = m_s[i+2] - s;

    Point A1_or_A3;
    Real denom = 1/(m_s[i] - m_s[i-1]);
    for(size_t j = 0; j < dimension; ++j)
    {
        A1_or_A3[j] = denom*(s1s*m_pnts[i-1][j] - s0s*m_pnts[i][j]);
    }

    Point A2_or_B2;
    for(size_t j = 0; j < dimension; ++j)
    {
        A2_or_B2[j] = denom21*(s2s*m_pnts[i][j] - s1s*m_pnts[i+1][j]);
    }

    Point B1_or_C;
    denom = 1/(m_s[i+1] - m_s[i-1]);
    for(size_t j = 0; j < dimension; ++j)
    {
        B1_or_C[j] = denom*(s2s*A1_or_A3[j] - s0s*A2_or_B2[j]);
    }

    denom = 1/(m_s[i+2] - m_s[i+1]);
    for(size_t j = 0; j < dimension; ++j)
    {
        A1_or_A3[j] = denom*(s3s*m_pnts[i+1][j] - s2s*m_pnts[i+2][j]);
    }

    Point B2;
    denom = 1/(m_s[i+2] - m_s[i]);
    for(size_t j = 0; j < dimension; ++j)
    {
        B2[j] = denom*(s3s*A2_or_B2[j] - s1s*A1_or_A3[j]);
    }

    for(size_t j = 0; j < dimension; ++j)
    {
        B1_or_C[j] = denom21*(s2s*B1_or_C[j] - s1s*B2[j]);
    }

    return B1_or_C;
}

template<class Real, class Point, size_t dimension>
Point catmull_rom<Real, Point, dimension>::prime(Real s) const
{
    // https://math.stackexchange.com/questions/843595/how-can-i-calculate-the-derivative-of-a-catmull-rom-spline-with-nonuniform-param
    // http://denkovacs.com/2016/02/catmull-rom-spline-derivatives/
    if (s < 0 || s > m_max_s)
    {
        throw std::domain_error("Parameter outside bounds.\n");
    }
    auto it = std::upper_bound(m_s.begin(), m_s.end(), s);
    //Now *it >= s. We want the index such that m_s[i] <= s < m_s[i+1]:
    size_t i = std::distance(m_s.begin(), it - 1);
    // We'll keep the assert in here a while until we feel good that we've understood this algorithm.
    assert(m_s[i] <= s && s < m_s[i+1]);
    Point A1;
    Real denom = 1/(m_s[i] - m_s[i-1]);
    Real k1 = (m_s[i]-s)*denom;
    Real k2 = (s - m_s[i-1])*denom;
    for (size_t j = 0; j < dimension; ++j)
    {
        A1[j] = k1*m_pnts[i-1][j] + k2*m_pnts[i][j];
    }

    Point A1p;
    for (size_t j = 0; j < dimension; ++j)
    {
        A1p[j] = denom*(m_pnts[i][j] - m_pnts[i-1][j]);
    }

    Point A2;
    denom = 1/(m_s[i+1] - m_s[i]);
    k1 = (m_s[i+1]-s)*denom;
    k2 = (s - m_s[i])*denom;
    for (size_t j = 0; j < dimension; ++j)
    {
        A2[j] = k1*m_pnts[i][j] + k2*m_pnts[i+1][j];
    }

    Point A2p;
    for (size_t j = 0; j < dimension; ++j)
    {
        A2p[j] = denom*(m_pnts[i+1][j] - m_pnts[i][j]);
    }


    Point B1;
    for (size_t j = 0; j < dimension; ++j)
    {
        B1[j] = k1*A1[j] + k2*A2[j];
    }

    Point A3;
    denom = 1/(m_s[i+2] - m_s[i+1]);
    k1 = (m_s[i+2]-s)*denom;
    k2 = (s - m_s[i+1])*denom;
    for (size_t j = 0; j < dimension; ++j)
    {
        A3[j] = k1*m_pnts[i+1][j] + k2*m_pnts[i+2][j];
    }

    Point A3p;
    for (size_t j = 0; j < dimension; ++j)
    {
        A3p[j] = denom*(m_pnts[i+2][j] - m_pnts[i+1][j]);
    }

    Point B2;
    denom = 1/(m_s[i+2] - m_s[i]);
    k1 = (m_s[i+2]-s)*denom;
    k2 = (s - m_s[i])*denom;
    for (size_t j = 0; j < dimension; ++j)
    {
        B2[j] = k1*A2[j] + k2*A3[j];
    }

    Point B1p;
    denom = 1/(m_s[i+1] - m_s[i-1]);
    for (size_t j = 0; j < dimension; ++j)
    {
        B1p[j] = denom*(A2[j] - A1[j] + (m_s[i+1]- s)*A1p[j] + (s-m_s[i-1])*A2p[j]);
    }

    Point B2p;
    denom = 1/(m_s[i+2] - m_s[i]);
    for (size_t j = 0; j < dimension; ++j)
    {
        B2p[j] = denom*(A3[j] - A2[j] + (m_s[i+2] - s)*A2p[j] + (s - m_s[i])*A3p[j]);
    }

    Point Cp;
    denom = 1/(m_s[i+1] - m_s[i]);
    for (size_t j = 0; j < dimension; ++j)
    {
        Cp[j] = denom*(B2[j] - B1[j] + (m_s[i+1] - s)*B1p[j] + (s - m_s[i])*B2p[j]);
    }
    return Cp;
}


}}
#endif
