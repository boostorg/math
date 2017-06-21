// Copyright Nick Thompson, 2017
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_QUADRATURE_DETAIL_SINH_SINH_DETAIL_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_SINH_SINH_DETAIL_HPP

#include <cmath>
#include <string>
#include <vector>
#include <typeinfo>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/next.hpp>

namespace boost{ namespace math{ namespace quadrature { namespace detail{


// Returns the sinh-sinh quadrature of a function f over the entire real line

template<class Real, class Policy>
class sinh_sinh_detail
{
public:
    sinh_sinh_detail(Real tol, size_t max_refinements);

    template<class F>
    Real integrate(const F f, Real* error, Real* L1) const;

private:
    Real m_tol;

    std::vector<std::vector<Real>> m_abscissas;
    std::vector<std::vector<Real>> m_weights;
};

template<class Real, class Policy>
sinh_sinh_detail<Real, Policy>::sinh_sinh_detail(Real tol, size_t max_refinements)
{
    using std::log;
    using std::sqrt;
    using std::cosh;
    using std::sinh;
    using std::ceil;
    using boost::math::constants::two_div_pi;
    using boost::math::constants::half_pi;

    m_tol = tol;

    // t_max is chosen to make g'(t_max) ~ sqrt(max) (g' grows faster than g).
    // This will allow some flexibility on the users part; they can at least square a number function without overflow.
    // But there is no unique choice; the further out we can evaluate the function, the better we can do on slowly decaying integrands.
    const Real t_max = log(2*two_div_pi<Real>()*log(2*two_div_pi<Real>()*sqrt(tools::max_value<Real>())));
    m_abscissas.resize(max_refinements);
    m_weights.resize(max_refinements);

    for (size_t i = 0; i < max_refinements; ++i)
    {
        Real h = (Real) 1/ (Real) (1<<i);
        size_t k = (size_t) ceil(t_max/(2*h));
        m_abscissas[i].reserve(k);
        m_weights[i].reserve(k);
        // We don't add 0 to the abscissas; that's treated as a special case.
        Real arg = h;
        while(arg < t_max)
        {
            Real tmp = half_pi<Real>()*sinh(arg);
            Real x = sinh(tmp);
            m_abscissas[i].emplace_back(x);
            Real w = cosh(arg)*half_pi<Real>()*cosh(tmp);
            m_weights[i].emplace_back(w);

            if (i != 0)
            {
                arg += 2*h;
            }
            else
            {
                arg += h;
            }
        }
    }
}

template<class Real, class Policy>
template<class F>
Real sinh_sinh_detail<Real, Policy>::integrate(const F f, Real* error, Real* L1) const
{
    using std::abs;
    using std::sqrt;
    using boost::math::constants::half;
    using boost::math::constants::half_pi;

    static const char* function = "boost::math::quadrature::sinh_sinh<%1%>::integrate";

    Real y_max = f(boost::math::tools::max_value<Real>());
    if(abs(y_max) > boost::math::tools::epsilon<Real>() || !(boost::math::isfinite)(y_max))
    {
        return policies::raise_domain_error(function,
           "The function you are trying to integrate does not go to zero at infinity, and instead evaluates to %1%", y_max, Policy());
    }

    Real y_min = f(std::numeric_limits<Real>::lowest());
    if(abs(y_min) > boost::math::tools::epsilon<Real>() || !(boost::math::isfinite)(y_min))
    {
        return policies::raise_domain_error(function,
           "The function you are trying to integrate does not go to zero at -infinity, and instead evaluates to %1%", y_max, Policy());
    }

    // Get the party started with two estimates of the integral:
    Real I0 = f(0)*half_pi<Real>();
    Real L1_I0 = abs(I0);
    for(size_t i = 0; i < m_abscissas[0].size(); ++i)
    {
        Real x = m_abscissas[0][i];
        Real yp = f(x);
        Real ym = f(-x);
        I0 += (yp + ym)*m_weights[0][i];
        L1_I0 += (abs(yp)+abs(ym))*m_weights[0][i];
    }

    // Uncomment the estimates to work the convergence on the command line.
    // std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    // std::cout << "First estimate : " << I0 << std::endl;
    Real I1 = I0;
    Real L1_I1 = L1_I0;
    for (size_t i = 0; i < m_abscissas[1].size(); ++i)
    {
        Real x= m_abscissas[1][i];
        Real yp = f(x);
        Real ym = f(-x);
        I1 += (yp + ym)*m_weights[1][i];
        L1_I1 += (abs(yp) + abs(ym))*m_weights[1][i];
    }

    I1 *= half<Real>();
    L1_I1 *= half<Real>();
    Real err = abs(I0 - I1);
    // std::cout << "Second estimate: " << I1 << " Error estimate at level " << 1 << " = " << err << std::endl;

    for(size_t i = 2; i < m_abscissas.size(); ++i)
    {
        I0 = I1;
        L1_I0 = L1_I1;

        I1 = half<Real>()*I0;
        L1_I1 = half<Real>()*L1_I0;
        Real h = (Real) 1/ (Real) (1 << i);
        Real sum = 0;
        Real absum = 0;

        Real abterm1 = 1;
        Real eps = std::numeric_limits<Real>::epsilon()*L1_I1;
        for(size_t j = 0; j < m_weights[i].size(); ++j)
        {
            Real x = m_abscissas[i][j];
            Real yp = f(x);
            Real ym = f(-x);
            sum += (yp + ym)*m_weights[i][j];
            Real abterm0 = (abs(yp) + abs(ym))*m_weights[i][j];
            absum += abterm0;

            // We require two consecutive terms to be < eps in case we hit a zero of f.
            if (x > (Real) 100 && abterm0 < eps && abterm1 < eps)
            {
                break;
            }
            abterm1 = abterm0;
        }

        I1 += sum*h;
        L1_I1 += absum*h;
        err = abs(I0 - I1);
        // std::cout << "Estimate:        " << I1 << " Error estimate at level " << i  << " = " << err << std::endl;
        if (!(boost::math::isfinite)(I1))
        {
            const char* err_msg = "The sinh_sinh quadrature evaluated your function at a singular point, leading to the value %1%.\n"
               "sinh_sinh quadrature cannot handle singularities in the domain.\n"
               "If you are sure your function has no singularities, please submit a bug against boost.math\n";
            return policies::raise_evaluation_error(function, err_msg, I1, Policy());
        }
        if (err <= m_tol*L1_I1)
        {
            break;
        }
    }

    if (error)
    {
        *error = err;
    }

    if (L1)
    {
        *L1 = L1_I1;
    }

    return I1;
}

}}}}
#endif
