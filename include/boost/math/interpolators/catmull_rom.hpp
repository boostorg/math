// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// This computes the Catmull-Rom spline from a list of points.

#ifndef BOOST_MATH_INTERPOLATORS_CATMULL_ROM
#define BOOST_MATH_INTERPOLATORS_CATMULL_ROM

#include <cmath>
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
    catmull_rom(const Point* const points, size_t num_pnts, bool closed = false, Real alpha = (Real) 1/ (Real) 2)
    {
        if (num_pnts < 4)
        {
            throw std::domain_error("The Catmull-Rom curve requires at least 4 points.\n");
        }
        using std::abs;
        m_s.resize(num_pnts+3);
        m_pnts.resize(num_pnts+3);
        if (closed)
        {
            m_pnts[0] = points[num_pnts-1];
            for (size_t i = 0; i < num_pnts; ++i) {
              m_pnts[i+1] = points[i];
            }
            m_pnts[num_pnts+1] = points[0];
            m_pnts[num_pnts+2] = points[1];
            m_s[0] = -detail::alpha_distance<Real, Point, dimension>(m_pnts[0], m_pnts[1], alpha);
            if (abs(m_s[0]) < std::numeric_limits<Real>::epsilon())
            {
                throw std::domain_error("The first and last point should not be the same for a closed curve-that is indicated by a boolean flag.\n");
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
            m_max_s = m_s[num_pnts];
        }
        else
        {
            for(size_t i = 0; i < dimension; ++i)
            {
                m_pnts[0][i] = 0;
                m_pnts[num_pnts+1][i] = 0;
                m_pnts[num_pnts+2][i] = 0;
            }

            for (size_t i = 0; i < num_pnts; ++i)
            {
                m_pnts[i+1] = points[i];
            }
            m_pnts[num_pnts+1] = points[0];
            m_pnts[num_pnts+2] = points[1];
            m_s[0] = -1;
            m_s[1] = 0;
            for (size_t i = 2; i < num_pnts+1; ++i)
            {
                Real d = detail::alpha_distance<Real, Point, dimension>(m_pnts[i-1], m_pnts[i-2], alpha);
                if (abs(d) < std::numeric_limits<Real>::epsilon())
                {
                    throw std::domain_error("The control points of the Catmull-Rom curve are too close together; this will lead to ill-conditioning.\n");
                }
                m_s[i] = m_s[i-1] + d;
            }
            m_max_s = m_s[num_pnts];
            m_s[num_pnts+2] = m_s[num_pnts+1] + 1;
        }
    }

    Real max_parameter() const
    {
        return m_max_s;
    }

    Point operator()(Real s) const
    {
        if (s < 0 || s > m_max_s)
        {
            throw std::domain_error("Parameter outside bounds.\n");
        }
        auto it = std::lower_bound(m_s.begin(), m_s.end(), s);
        size_t i = std::distance(m_s.begin(), it);
        //std::cout << "Index = " << i << std::endl;
        assert(m_s[i] <= s && s < m_s[i+1]);
        Point A1;
        Real denom = 1/(m_s[i] - m_s[i-1]);
        Real k1 = (m_s[i]-s)*denom;
        Real k2 = (s - m_s[i-1])*denom;
        for (size_t j = 0; j < dimension; ++j)
        {
            A1[j] = k1*m_pnts[i-1][j] + k2*m_pnts[i][j];
        }

        Point A2;
        denom = 1/(m_s[i+1] - m_s[i]);
        k1 = (m_s[i+1]-s)*denom;
        k2 = (s - m_s[i])*denom;
        for (size_t j = 0; j < dimension; ++j)
        {
            A2[j] = k1*m_pnts[i][j] + k2*m_pnts[i+1][j];
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

        Point B2;
        denom = 1/(m_s[i+2] - m_s[i]);
        k1 = (m_s[i+2]-s)*denom;
        k2 = (s - m_s[i])*denom;
        for (size_t j = 0; j < dimension; ++j)
        {
            B2[j] = k1*A2[j] + k2*A3[j];
        }

        Point C;
        denom = 1/(m_s[i+1] - m_s[i]);
        k1 = (m_s[i+1]-s)*denom;
        k2 = (s - m_s[i])*denom;
        for (size_t j = 0; j < dimension; ++j)
        {
            C[j] = k1*B1[j] + k2*B2[j];
        }


        return C;
    }

private:
    std::vector<Point> m_pnts;
    std::vector<Real> m_s;
    Real m_max_s;
};

}}
#endif
