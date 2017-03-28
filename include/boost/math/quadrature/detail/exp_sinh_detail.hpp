#ifndef BOOST_MATH_QUADRATURE_DETAIL_EXP_SINH_DETAIL_HPP
#define BOOST_MATH_QUADRATURE_DETAIL_EXP_SINH_DETAIL_HPP

#include <cmath>
#include <vector>
#include <typeinfo>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/next.hpp>
#include <boost/format.hpp>

namespace boost{ namespace math{ namespace detail{


// Returns the exp-sinh quadrature of a function f over the open interval (0, infinity)

template<class Real>
class exp_sinh_detail
{
public:
    exp_sinh_detail(Real tol, size_t max_refinements);

    template<class F>
    Real integrate(const F f, Real* error, Real* L1) const;

private:
    Real m_tol;

    std::vector<std::vector<Real>> m_abscissas;
    std::vector<std::vector<Real>> m_weights;
};

template<class Real>
exp_sinh_detail<Real>::exp_sinh_detail(Real tol, size_t max_refinements)
{
    using std::exp;
    using std::log;
    using std::sqrt;
    using std::asinh;
    using std::cosh;
    using std::sinh;
    using std::ceil;
    using boost::math::constants::two_div_pi;
    using boost::math::constants::half_pi;

    m_tol = tol;
    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);
    // t_min is chosen such that x = exp(pi/2 sinh(-t_max)) = eps/10000.
    // This is a compromise; using some denormals so as to pick up information about singularities at zero,
    // while not risking hitting the singularity through roundoff error.
    const Real t_min = asinh(two_div_pi<Real>()*log(std::numeric_limits<Real>::epsilon()/10000));

    // t_max is chosen to make g'(t_max) ~ sqrt(max) (g' grows faster than g).
    // This will allow some flexibility on the users part; they can at least square a number function without overflow.
    // But there is no unique choice; the further out we can evaluate the function, the better we can do on slowly decaying integrands.
    const Real t_max = log(2*two_div_pi<Real>()*log(2*two_div_pi<Real>()*sqrt(std::numeric_limits<Real>::max())));
    m_abscissas.resize(max_refinements);
    m_weights.resize(max_refinements);

    for (size_t i = 0; i < max_refinements; ++i)
    {
        Real h = (Real) 1/ (Real) (1<<i);
        size_t k = (size_t) ceil((t_max - t_min)/(2*h));
        m_abscissas[i].reserve(k);
        m_weights[i].reserve(k);
        Real arg = t_min;
        if (i != 0)
        {
            arg = t_min + h;
        }

        while(arg < t_max)
        {
            Real x = exp(half_pi<Real>()*sinh(arg));
            m_abscissas[i].emplace_back(x);
            Real w = cosh(arg)*half_pi<Real>()*x;
            m_weights[i].emplace_back(w);

            if (i != 0)
            {
                arg = arg + 2*h;
            }
            else
            {
                arg = arg + h;
            }
        }
    }
}

template<class Real>
template<class F>
Real exp_sinh_detail<Real>::integrate(const F f, Real* error, Real* L1) const
{
    using std::abs;
    using std::floor;
    using std::tanh;
    using std::sinh;
    using std::sqrt;
    using boost::math::constants::half;
    using boost::math::constants::half_pi;

    Real y_max = f(std::numeric_limits<Real>::max());
    if(abs(y_max) > std::numeric_limits<Real>::epsilon() || !isfinite(y_max))
    {
        boost::basic_format<char> err_msg = boost::format("\nThe function you are trying to integrate does not go to zero at infinity, and instead evaluates to %1% at x=%2%.\n") % y_max % std::numeric_limits<Real>::max();
        throw std::domain_error(err_msg.str());
    }

    std::cout << std::setprecision(std::numeric_limits<Real>::digits10);

    // Get the party started with two estimates of the integral:
    Real I0 = 0;
    Real L1_I0 = 0;
    for(size_t i = 0; i < m_abscissas[0].size(); ++i)
    {
        Real y = f(m_abscissas[0][i]);
        I0 += y*m_weights[0][i];
        L1_I0 += abs(y)*m_weights[0][i];
    }

    // std::cout << "First estimate : " << I0 << std::endl;
    Real I1 = I0;
    Real L1_I1 = L1_I0;
    for (size_t i = 0; i < m_abscissas[1].size(); ++i)
    {
        Real y = f(m_abscissas[1][i]);
        I1 += y*m_weights[1][i];
        L1_I1 += abs(y)*m_weights[1][i];
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
            Real y = f(x);
            sum += y*m_weights[i][j];
            Real abterm0 = abs(y)*m_weights[i][j];
            absum += abterm0;

            // We require two consecutive terms to be < eps in case we hit a zero of f.
            // Numerical experiments indicate that we should start this check at ~30 for floats,
            // ~60 for doubles, and ~100 for long doubles.
            // However, starting the check at x = 10 rather than x = 100 will only save two function evaluations.
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
        if (!isfinite(I1))
        {
            throw std::domain_error("The exp_sinh quadrature evaluated your function at a singular point. Please ensure your function evaluates to a finite number of its entire domain.\n");
        }
        if (err <= m_tol*L1_I1 || err <= std::numeric_limits<Real>::epsilon() /* catch the 0 case. */)
        {
            break;
        }
    }

    if (error)
    {
        *error = err;
    }

    if(L1)
    {
        *L1 = L1_I1;
    }

    if(L1_I1 != (Real) 0 && L1_I1/abs(I1) > 1/std::numeric_limits<Real>::epsilon())
    {
        Real cond = abs(I1)/L1_I1;
        Real inv_prec = 1.0/std::numeric_limits<Real>::epsilon();
        boost::basic_format<char> err_msg = boost::format("\nThe condition number of the quadrature sum is %1%, which exceeds the inverse of the precision of the type %2%.\nNo correct digits can be expected for the integral.\n") % cond % inv_prec;
        throw boost::math::evaluation_error(err_msg.str());
    }

    return I1;
}

}}}
#endif
