// Copyright Nick Thompson, 2019
// Use, modification and distribution are subject to the
// Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

/*
 * References:
 * Ooura, Takuya. "A double exponential formula for the Fourier transforms." Publications of the Research Institute for Mathematical Sciences 41.4 (2005): 971-977.
 * Ooura, Takuya, and Masatake Mori. "A robust double exponential formula for Fourier-type integrals." Journal of computational and applied mathematics 112.1-2 (1999): 229-241.
 * http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html
 */
#ifndef BOOST_MATH_QUADRATURE_OOURA_FOURIER_INTEGRALS_HPP
#define BOOST_MATH_QUADRATURE_OOURA_FOURIER_INTEGRALS_HPP
#include <utility>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/constants/constants.hpp>

namespace boost { namespace math { namespace quadrature {

namespace detail {
    template<class Real>
    Real ooura_substitution(Real t, Real alpha)
    {
        using std::expm1;
        Real beta = (Real) 1/ (Real) 4;
        if (t == 0)
        {
            return 1/(2 + alpha + beta);
        }
        Real z1 = -beta*expm1(t);
        Real z2 = alpha*expm1(-t);
        Real theta = -2*t +z1 + z2;
        return -t/expm1(theta);
    }

    template<class Real>
    Real ooura_substitution_prime(Real t, Real alpha)
    {
        Real beta = (Real) 1/ (Real) 4;
        if(t==0)
        {
            Real num = 4 + alpha*alpha + 3*beta + beta*beta + alpha*(5+2*beta);
            Real denom = 2*(2+alpha+beta)*(2+alpha+beta);
            return num/denom;
        }
        using std::exp;
        using std::expm1;
        using std::abs;
        using std::pow;
        using std::sqrt;
        using std::numeric_limits;
        Real z1 = -beta*expm1(t);
        Real z2 = alpha*expm1(-t);
        Real theta = -2*t +z1 + z2;
        Real B = -expm1(theta);
        if (abs(B) > sqrt(numeric_limits<Real>::max()))
        {
            return 0;
        }
        Real inv_B = 1/B;
        Real et = exp(t);
        Real emt = exp(-t);
        Real dBdt = (2+beta*et+alpha*emt)*exp(theta);
        Real res = inv_B*(1-t*dBdt*inv_B);
        return inv_B*(1-t*dBdt*inv_B);
    }

    template<class Real>
    std::pair<Real, Real> ooura_pair(Real t, Real alpha)
    {
        using std::expm1;
        using std::exp;
        using std::abs;
        using std::pow;
        using std::sqrt;
        using std::numeric_limits;

        Real beta = (Real) 1/ (Real) 4;
        if (t == 0)
        {
            Real num = 4 + alpha*alpha + 3*beta + beta*beta + alpha*(5+2*beta);
            Real denom = 2*(2+alpha+beta)*(2+alpha+beta);
            return std::make_pair<Real, Real>(1/(2 + alpha + beta), num/denom);
        }
        Real etm1 = expm1(t);
        Real z1 = -beta*etm1;
        Real emtm1 = expm1(-t);
        Real z2 = alpha*emtm1;
        Real theta = -2*t +z1 + z2;
        Real B = -expm1(theta);
        Real phi = t/B;

        if (abs(B) > sqrt(numeric_limits<Real>::max()))
        {
            return {phi, 0};
        }
        Real inv_B = 1/B;
        Real dBdt = (2+beta*(etm1+1)+alpha*(emtm1+1))*(-B+1);
        Real res = inv_B*(1-t*dBdt*inv_B);
        Real phi_prime =  inv_B*(1-t*dBdt*inv_B);
        return {phi, phi_prime};
    }


    template<class Real>
    Real determine_alpha(Real M)
    {
        using boost::math::constants::pi;
        using std::log1p;
        using std::sqrt;
        Real x = sqrt(1 + M*log1p(M)/(4*pi<Real>()) );
        return 1/(4*x);
    }

    template<class Real>
    Real find_max_sensible_t(Real alpha)
    {
        using std::nextafter;
        using std::numeric_limits;
        Real eps = numeric_limits<Real>::epsilon();
        Real t = 3;
        Real phi;
        Real diff = 0.5;

        while (diff > t*numeric_limits<Real>::epsilon())
        {
            do
            {
                t += diff;
                phi = ooura_substitution<Real>(t, alpha);
            } while(phi != t);
            t = t - diff;
            diff /= 2;
        }
        return t;
    }


    template<class Real>
    Real find_min_sensible_t(Real alpha)
    {
        using std::nextafter;
        using std::numeric_limits;
        using std::min;
        using std::abs;
        Real eps = numeric_limits<Real>::epsilon();
        Real t = -1;
        Real phi;
        Real phi_prime;
        Real x;
        Real diff = 0.5;
        while (diff > abs(t*std::numeric_limits<Real>::epsilon()))
        {
            do
            {
                t -= diff;
                phi = ooura_substitution<Real>(t, alpha);
                phi_prime = ooura_substitution_prime<Real>(t, alpha);
                x = min(phi, phi_prime);
            } while(x != 0);
            t += diff;
            diff /= 2;
        }
        return t;
    }
}

template<class F, class Real>
Real ooura_fourier_sin(F const & integrand, const Real omega, const Real tolerance = tools::root_epsilon<Real>())
{
    using std::max;
    using std::min;
    using std::sin;
    using std::abs;
    using std::numeric_limits;
    using boost::math::constants::pi;
    if (omega == 0)
    {
        return 0;
    }
    if (omega < 0)
    {
        return -ooura_fourier_sin(integrand, -omega, tolerance);
    }

    Real h = Real(1)/Real(2);
    Real I0 = 0;
    Real I1 = numeric_limits<Real>::max();
    Real error_estimate = numeric_limits<Real>::quiet_NaN();
    int k = 0;
    do {
        Real M = pi<Real>()/h;
        Real alpha = detail::determine_alpha(M);

        Real t;
        Real phi;
        long long n = 0;
        int consecutive = 0;
        do {
            t = n*h;
            auto p = detail::ooura_pair(t, alpha);
            phi = p.first;
            Real b = p.second*sin(M*phi);
            Real arg = M*phi/omega;
            Real term;
            if (arg == 0)
            {
                // then b = 0:
                term = 0;
            }
            else
            {
                term = integrand(M*phi/omega)*b;
            }
            if (I0 != 0 && abs(term/I0) < numeric_limits<Real>::epsilon())
            {
                ++consecutive;
            }
            else
            {
                consecutive = 0;
            }
            I0 += term;
            ++n;
        } while(phi != t && consecutive < 10);

        n = 1;
        consecutive = 0;
        do {
            t = -n*h;
            auto p = detail::ooura_pair(t, alpha);
            phi = p.first;
            Real b = p.second*sin(M*phi);
            Real arg = M*phi/omega;
            Real term;
            if (arg == 0)
            {
                term = 0;
            }
            else
            {
                term = integrand(arg)*b;
            }
            if (I0 != 0 && abs(term/I0) < numeric_limits<Real>::epsilon())
            {
                ++consecutive;
            }
            else
            {
                consecutive = 0;
            }
            I0 += term;
            ++n;
        } while(phi > numeric_limits<Real>::min() && consecutive < 10);

        h /= 4;
        error_estimate = abs(I0 - I1)/max( abs(I0), abs(I1));
        I1 = I0;
        I0 = 0;
    } while (error_estimate > tolerance && ++k < 10);

    return pi<Real>()*I1/omega;
}


template<class F, class Real>
Real ooura_fourier_cos(F const & integrand, const Real omega, const Real tolerance = tools::root_epsilon<Real>())
{
    using std::max;
    using std::min;
    using std::sin;
    using std::abs;
    using std::numeric_limits;
    using boost::math::constants::pi;
    using boost::math::constants::half_pi;
    if (omega == 0)
    {
        throw std::domain_error("At omega = 0, the integral is not oscillatory. The user must choose an appropriate method for this case.\n");
    }
    if (omega < 0)
    {
        return ooura_fourier_cos(integrand, -omega, tolerance);
    }

    Real h = Real(1)/Real(2);
    Real I0 = 0;
    Real I1 = numeric_limits<Real>::max();
    Real error_estimate = numeric_limits<Real>::quiet_NaN();
    do {
        Real M = pi<Real>()/h;
        Real alpha = detail::determine_alpha(M);

        Real t;
        Real phi;
        long long n = 0;
        int consecutive = 0;
        do {
            t = n*h - half_pi<Real>()/M;
            auto p = detail::ooura_pair(t, alpha);
            phi = p.first;
            Real b = p.second*cos(M*phi);
            Real arg = M*phi/omega;
            Real term;
            if (arg == 0)
            {
                // then b = 0:
                term = 0;
            }
            else
            {
                term = integrand(M*phi/omega)*b;
            }
            if (I0 != 0 && abs(term/I0) < numeric_limits<Real>::epsilon())
            {
                ++consecutive;
            }
            else
            {
                consecutive = 0;
            }
            I0 += term;
            ++n;
        } while(phi != t && consecutive < 10);

        n = 1;
        consecutive = 0;
        do {
            t = -n*h - half_pi<Real>()/M;
            auto p = detail::ooura_pair(t, alpha);
            phi = p.first;
            Real b = p.second*cos(M*phi);
            Real arg = M*phi/omega;
            Real term;
            if (arg == 0)
            {
                term = 0;
            }
            else
            {
                term = integrand(arg)*b;
            }
            if (I0 != 0 && abs(term/I0) < numeric_limits<Real>::epsilon())
            {
                ++consecutive;
            }
            else
            {
                consecutive = 0;
            }
            I0 += term;
            ++n;
        } while(phi > numeric_limits<Real>::min() && consecutive < 10);

        h /= 2;
        error_estimate = abs(I0 - I1)/max( abs(I0), abs(I1));
        I1 = I0;
        I0 = 0;
    } while (error_estimate > tolerance);

    return pi<Real>()*I1/omega;
}


}}}
#endif
