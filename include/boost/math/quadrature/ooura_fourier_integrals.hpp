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
#include <memory>
#include <boost/math/quadrature/detail/ooura_fourier_integrals_detail.hpp>

namespace boost { namespace math { namespace quadrature {

template<class Real>
class ooura_fourier_sin {
public:
    ooura_fourier_sin(const Real relative_error_tolerance = tools::root_epsilon<Real>(), size_t levels = sizeof(Real)) : impl_(std::make_shared<detail::ooura_fourier_sin_detail<Real>>(relative_error_tolerance, levels))
    {}

    template<class F>
    Real integrate(F const & f, Real omega) {
        return impl_->integrate(f, omega);
    }
private:
    std::shared_ptr<detail::ooura_fourier_sin_detail<Real>> impl_;
};

/*
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
}*/


}}}
#endif
